# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
import lingpy
from lingpy.util import PROG
import argparse
from six import text_type as str
from six import with_metaclass

class CommandMeta(type):
    """
    A metaclass which keeps track of subclasses, if they have all-lowercase names.
    """
    __instances = []

    def __init__(self, name, bases, dct):
        super(CommandMeta, self).__init__(name, bases, dct)
        if name == name.lower():
            self.__instances.append(self)

    def __iter__(self):
        return iter(self.__instances)


class Command(with_metaclass(CommandMeta, object)):
    """Base class for subcommands of the lingpy command line interface."""
    help = None

    @classmethod
    def subparser(cls, parser):
        """Hook to define subcommand arguments."""
        return

    def output(self, args, content):
        if args.output_file:
            lingpy.util.write_text_file(args.output_file, content)
        else:
            print(content)

    def __call__(self, args):
        """Hook to run the subcommand."""
        raise NotImplemented()


def _cmd_by_name(name):
    for cmd in Command:
        if cmd.__name__ == name:
            return cmd()


def add_option(parser, name_, default_, help_, short_opt=None, **kw):
    names = ['--' + name_]
    if short_opt:
        names.append('-' + short_opt)

    if isinstance(default_, bool):
        kw['action'] = 'store_false' if default_ else 'store_true'
    elif isinstance(default_, int):
        kw['type'] = float
    elif isinstance(default_, float):
        kw['type'] = float
    kw['default'] = default_
    kw['help'] = help_
    parser.add_argument(*names, **kw)


def add_shared_args(p):
    add_option(p, 'input-file', None, "Path to input file.", short_opt='i')
    add_option(p, 'output-file', None, "Path to output file.", short_opt='o')
    add_option(p, 'factor', 0.3, "Factor for prosodic strings in SCA alignments.")
    add_option(p, 'gop', -1.0, "Gap opening penalty.")
    add_option(p, 'scale', 0.5, "Gap extension scale.")
    add_option(
        p,
        'restricted-chars',
        'T_',
        "Characters in the sound classes which mark word or morpheme boundaries.")


def add_cognate_identifier_option(p, default):
    add_option(
        p,
        'cognate-identifier',
        default,
        'Name for the column with the cognate judgments.',
        short_opt='c')


def add_tree_calc_option(p):
    add_option(
        p,
        'tree-calc',
        'upgma',
        "Select the tree cluster method you want to use for the guide tree.",
        choices=['ugpma', 'neighbor'])


def add_method_option(p, default, choices, spec=''):
    add_option(
        p,
        'method',
        default,
        "The %s method you want to use." % spec,
        choices=choices)

def add_format_option(p, default, choices):
    add_option(p, 'format', default, "Output format.", choices=choices)


def add_strings_option(p, n):
    add_option(
        p, 'strings', None, "Input sequences you want to align.", nargs=n, short_opt='s')


def add_mode_option(p, choices):
    add_option(p, 'mode', 'global', 'Alignment mode', choices=choices, short_opt='m')


def add_align_method_option(p):
    add_option(
        p,
        'align-method',
        'library',
        "Select how you want to align, based on a library (T-COFFEE), or " +
        "in a classical progressive way.",
        choices=['library', 'progressive'])


class wordlist(Command):
    """
    Load a wordlist and carry out simple checks.
    """
    @classmethod
    def subparser(cls, p):
        add_shared_args(p)
        add_option(p, 'check', False, 'Check the content of your wordlist.')
        add_option(p, 'stats', False, 'Retrieve basic statistics of your wordlist.')
        add_option(p, 'transform', False, 'Transform your wordlist.', short_opt='t')
        add_cognate_identifier_option(p, 'cogid')
        add_option(
            p,
            'cluster-method',
            'upgma',
            'Cluster method you want to use.',
            choices=['upgma', 'single', 'complete', 'mcl', 'ward', 'link_clustering'])
        add_option(
            p,
            'columns',
            None,
            "Select the columns and the order in which you want to write your wordlist.",
            nargs='+')
        add_option(
            p,
            'add-row',
            ['ipa', 'tokens', 'ipa2tokens'],
            'Add rows to wordlist by converting source with target via function.',
            nargs=3)
        add_option(
            p,
            'calculate',
            '',
            "Calculate a tree from the wordlist.",
            choices=['dst', 'tree', 'diversity'])
        add_format_option(
            p,
            'tsv',
            ['dst', 'taxa', 'paps.nex', 'nwk', 'cluster', 'starling', 'multistate.nex',
             'separated'])
        add_option(
            p,
            'tree-calc',
            'upgma',
            "Select the tree cluster method you want to use for the guide tree.",
            choices=['ugpma', 'neighbor'])
        add_option(p, 'missing', '?', "Missing value for the output to Nexus files.")

    def __call__(self, args):
        if args.input_file:
            if args.check:
                wl = lingpy.compare.lexstat.LexStat(args.input_file, check=True)
            else:
                wl = lingpy.basic.wordlist.Wordlist(args.input_file)

            # process stuff according to arguments, will surely need to be further
            # refined
            if args.stats:
                print('---BASIC---')
                print('Height:  {0}'.format(wl.height))
                print('Width:   {0}'.format(wl.width))
                print('Length:  {0}'.format(len(wl)))
                print('')
                print('---COVERAGE---')
                for k, v in wl.coverage(stats='absolute').items():
                    print('{0:20}  :  {1}'.format(k, v))
                print('')
            if args.calculate:
                wl.calculate(args.calculate, ref=args.cognate_identifier,
                             tree_calc=args.tree_calc)
                if not args.output_file:
                    if args.calculate == 'tree':
                        print('---TREE---')
                        print(wl.tree.asciiArt())
                    if args.calculate == 'diversity':
                        print('---DIVERSITY---')
                        print(wl.diversity)

            if args.output_file:
                wl.output(
                    args.format, filename=args.output_file,
                    ref=args.cognate_identifier, missing=args.missing)


class alignments(Command):
    """
    Carry out alignment analysis of a wordlist file with readily detected cognates.
    """
    @classmethod
    def subparser(cls, p):
        add_shared_args(p)
        add_option(
            p,
            'mode',
            'global',
            "Select the alignment mode.",
            choices=['global', 'overlap', 'dialign'],
            short_opt='m')
        add_method_option(p, 'basic', ['sca', 'basic'])
        add_tree_calc_option(p)
        add_option(p, 'gap-weight', 0.5, 'Gap weight.')
        add_format_option(p, 'tsv', ['msa', 'html', 'tsv'])
        add_align_method_option(p)
        add_option(p, 'iterate', False, "Postprocess using iteration methods.")
        add_cognate_identifier_option(p, 'lingpyid')
        add_option(
            p,
            'alignment-identifier',
            'alignment',
            'Name for the column with the alignments.')
        add_option(
            p,
            'use-logodds',
            False,
            'Use logodds from LexStat to calculate the alignments.')

    def __call__(self, args):
        if args.input_file:
            alms = lingpy.align.sca.Alignments(
                args.input_file,
                merge_vowels=not args.no_merged_vowels,
                expand_nasals=args.expand_nasals,
                ref=args.cognate_identifier,
                alignment=args.alignment_identifier)
            if args.use_logodds:
                lex = lingpy.compare.lexstat.LexStat(args.input_file)
                if not hasattr(lex, 'cscorer'):
                    raise ValueError("No scorer has been submitted along with your file.")
            alms.align(method=args.align_method, gop=args.gop, scale=args.scale,
                       factor=args.factor, iteration=args.iterate, model=args.model,
                       tree_calc=args.tree_calc,
                       restricted_chars=args.restricted_chars, mode=args.mode)
            alms.output(args.format, filename=args.output_file, ignore='all')
        else:
            print('an input file is required')


class settings(Command):
    @classmethod
    def subparser(cls, p):
        add_option(
            p,
            'params',
            None,
            "Parameter names to show default settings for.",
            nargs='+',
            short_opt='p')

    def __call__(self, args):
        for k, v in lingpy.settings.rcParams.items():
            if not args.params or k in args.params:
                print('{0:20} : {1}'.format(k, repr(v)))


class lexstat(Command):
    @classmethod
    def subparser(cls, p):
        add_shared_args(p)
        add_method_option(
            p,
            'sca',
            ['lexstat', 'sca', 'edit-dist', 'turchin'],
            spec='cognate detection')
        add_option(p, 'runs', 1000, 'Number of iterations for the permutation test.')
        add_option(p, 'check', False, 'Carry out a check before running your analysis.')
        add_option(
            p,
            'scoring-threshold',
            lingpy.settings.rcParams['lexstat_scoring_threshold'],
            'Select the threshold to be used when creating the scorer.')
        add_option(
            p,
            'threshold',
            0.45,
            'Select the threshold to be used when creating the scorer.',
            short_opt='t')
        add_cognate_identifier_option(p, 'lingpyid')
        add_option(
            p,
            'cluster-method',
            'upgma',
            'Select the name of the cluster method you want to use.',
            choices=['upgma', 'single', 'complete', 'mcl', 'ward', 'link_clustering'])
        add_option(
            p,
            'ratio',
            [2, 1],
            'Ratio between the logodds and the sound class scorer.',
            nargs=2,
            type=int)

    def __call__(self, args):
        if args.input_file:
            lex = lingpy.compare.lexstat.LexStat(
                args.input_file,
                merge_vowels=not args.no_merged_vowels,
                expand_nasals=args.expand_nasals,
                check=args.check
            )
            if args.verbose:
                print("[i] Loaded file {0}.".format(args.input_file))

            if args.method == 'lexstat':
                # check for existing scorer
                if not hasattr(lex, 'cscorer'):
                    lex.get_scorer(runs=args.runs, preprocessing=False,
                                   ratio=args.ratio, factor=args.factor,
                                   threshold=args.scoring_threshold,
                                   restricted_chars=args.restricted_chars)
                    if args.verbose:
                        print("[i] Calculated the scorer.")
                else:
                    if args.verbose:
                        print(
                            "[i] Scorer has already been calculated. Skipping re-calculation.")
            lex.cluster(args.method, threshold=args.threshold,
                        cluster_method=args.cluster_method, scale=args.scale,
                        factor=args.factor, restricted_chars=args.restricted_chars,
                        gop=args.gop, ref=args.cognate_identifier)
            if args.verbose:
                print("[i] Clustered the data.")
            lex.output('tsv', filename=args.output_file)
            return len(lex.get_etymdict(ref=args.cognate_identifier))
        else:
            print('an input file is required')


class pairwise(Command):
    """
    Run pairwise analyses from command line in LingPy

    Notes
    -----

    Currently, the following options are supported:

    * run normal analyses without sound class strings
    * run sound-class based analyses

    Furthermore, input output is handled as follows:

    * define user input using psa-formats in lingpy
    * define user output (stdout, file)
    """
    @classmethod
    def subparser(cls, p):
        add_shared_args(p)
        add_strings_option(p, 2)
        add_mode_option(p, ['global', 'local', 'overlap', 'dialign'])
        add_option(
            p,
            'distance',
            False,
            "Choose whether you want distances or similarities to be reported.",
            short_opt='d')
        add_method_option(p, 'basic', ['sca', 'basic'], spec='basic')

    def __call__(self, args):
        def make_out(x, y, z):
            try:
                return '\t'.join(x) + '\n' + '\t'.join(y) + '\n{0:.2f}'.format(z)
            except TypeError:
                out1 = ''.join(x[0])+'\t|\t'+'\t'.join(x[1])+'\t|\t'+''.join(x[2])
                out2 = ''.join(y[0])+'\t|\t'+'\t'.join(y[1])+'\t|\t'+''.join(y[2])
                return out1+'\n'+out2+'\n'+'{0:.2f}'.format(z)

        if args.strings:
            if args.method == 'basic':
                almA, almB, sim = lingpy.align.pairwise.pw_align(
                    args.strings[0], args.strings[1], **args.__dict__)
            else:
                pair = lingpy.align.pairwise.Pairwise(
                    *args.strings,
                    merge_vowels=not args.no_merged_vowels,
                    expand_nasals=args.expand_nasals)
                pair.align(**args.__dict__)
                almA, almB, sim = pair.alignments[0]

            self.output(args, make_out(almA, almB, sim))
        elif args.input_file:
            pairs = lingpy.align.sca.PSA(args.input_file)
            output = ''
            if args.method == 'basic':
                for i, (seqA, seqB) in enumerate(pairs.tokens):
                    almA, almB, sim = lingpy.align.pairwise.pw_align(
                        ''.join(seqA), ''.join(seqB), **args.__dict__)
                    pairs.alignments[i] = [almA, almB, sim]
                    output += make_out(almA, almB, sim) + '\n'
            elif args.method == 'sca':
                pairs.align(**args.__dict__)
                for almA, almB, sim in pairs.alignments:
                    output += make_out(almA, almB, sim) + '\n'

            if args.output_file:
                pairs.output('psa', filename=args.output_file)
            else:
                print(output)
        else:
            print('either strings or an input file must be specified')


class multiple(Command):
    """
    Multiple alignment console interface for LingPy.

    Notes
    -----

    """
    @classmethod
    def subparser(cls, p):
        add_shared_args(p)
        add_strings_option(p, '+')
        add_mode_option(p, ['global', 'overlap', 'dialign'])
        add_method_option(p, 'basic', ['sca', 'basic'], spec='basic')
        add_tree_calc_option(p)
        add_option(
            p,
            'gap-weight',
            0.5,
            "Select the tree cluster method you want to use for the guide tree.")
        add_format_option(p, 'msa', ['msa', 'html', 'tex', 'psa'])
        add_align_method_option(p)
        add_option(
            p,
            'iteration',
            [],
            "Select method for postprocessing.",
            nargs="+",
            choices=['orphans', 'clusters', 'similar-gap-sites', 'all-sequences'])
        add_option(p, 'swap-check', False, "Search for swapped sites automatically.")

    def __call__(self, args):
        def align_msa(msa, args):
            if args.method == 'library':
                msa.lib_align(**args.__dict__)
            else:
                msa.prog_align(**args.__dict__)
            if 'orphans' in args.iteration:
                msa.iterate_orphans()
            if 'clusters' in args.iteration:
                msa.iterate_clusters()
            if 'similar-gap-sites' in args.iteration:
                msa.iterate_similar_gap_sites()
            if 'all-sequences' in args.iteration:
                msa.iterate_all_sequences()
            if args.swap_check:
                msa.swap_check()

        if args.strings:
            if args.method == 'basic':
                alms = lingpy.align.multiple.mult_align(
                    args.strings,
                    pprint=False,
                    gop=args.gop,
                    tree_calc=args.tree_calc,
                    scoredict=False,
                    scale=args.scale)
            elif args.method == 'sca':
                msa = lingpy.align.multiple.Multiple(
                    args.strings,
                    merge_vowels=not args.no_merged_vowels,
                    expand_nasals=args.expand_nasals,
                    model=args.model)
                align_msa(msa, args)
                alms = msa.alm_matrix

            self.output(args, '\n'.join(['\t'.join(seq) for seq in alms]))

        elif args.input_file:
            msa = lingpy.align.sca.MSA(
                args.input_file,
                merge_vowels=not args.no_merged_vowels,
                expand_nasals=args.expand_nasals,
                model=args.model)
            if args.method == 'basic':
                seqs = [str(''.join(t)) for t in msa.tokens]
                alms = lingpy.align.multiple.mult_align(
                    seqs,
                    pprint=False,
                    gop=args.gop,
                    tree_calc=args.tree_calc,
                    scoredict=False,
                    scale=args.scale,
                )
                msa.alm_matrix = []
                for alm in alms:
                    msa.alm_matrix += [alm]
            else:
                align_msa(msa, args)
            if args.output_file:
                msa.output(args.format, filename=args.output_file)
            else:
                for taxon, alm in zip(msa.taxa, msa.alm_matrix):
                    print('{0}\t{1}'.format(taxon, '\t'.join(alm)))
            return msa.taxa, msa.alm_matrix
        else:
            print('either strings or an input file must be specified')


class help(Command):
    """
    Show help for commands.
    """
    @classmethod
    def subparser(self, parser):
        parser.add_argument(
            'cmd', choices=[cmd.__name__ for cmd in Command if cmd.__name__ != 'help'])

    def __call__(self, args):
        cmd = _cmd_by_name(args.cmd)
        if cmd.__doc__:
            print('\n%s\n' % cmd.__doc__.strip())
        print(cmd.help)


def get_parser():
    # basic parser for lingpy
    parser = argparse.ArgumentParser(
        description=main.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_option(parser, 'version', None, '', action='version', version=PROG)
    add_option(
        parser,
        'schema',
        'ipa',
        "Input schema of your transcription system.",
        choices=['ipa', 'asjp'])
    add_option(
        parser,
        'model',
        'sca',
        'Select your sound class model.',
        choices=['sca', 'dolgo', 'asjp'])
    add_option(
        parser,
        'no-merged-vowels',
        False,
        'Do not merge vowels in automatic segmentation.')
    add_option(
        parser,
        'expand-nasals',
        False,
        "Display nasals in vowels as an extra segment in automatic segmentation.")
    add_option(parser, 'verbose', False, 'Trigger verbose output.', short_opt='v')

    subparsers = parser.add_subparsers(dest="subcommand")
    for cmd in Command:
        subparser = subparsers.add_parser(
            cmd.__name__,
            help=(cmd.__doc__ or '').strip().split('\n')[0],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        cmd.subparser(subparser)
        cmd.help = subparser.format_help()

    return parser


def main(*args):
    """
    LingPy command line interface.
    """
    args = get_parser().parse_args(args or None)
    lingpy.rc(schema=args.schema)
    lingpy.rc(model=lingpy.settings.rcParams[args.model])
    return _cmd_by_name(args.subcommand)(args)
