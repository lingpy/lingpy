# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
import lingpy
import argparse
from six import text_type as str

# basic parser for lingpy
parser = argparse.ArgumentParser(description='LingPy command line interface.')

parser.add_argument('--version', action='version', version=lingpy.__version__)
parser.add_argument('--schema', action='store', choices=['ipa', 'asjp'],    
        default='ipa', help="Select the input schema of your transcription system.")
parser.add_argument('--model', action='store', choices=['sca', 'dolgo', 'asjp'],
        default='sca', help='Select your sound class model.')
parser.add_argument('--merge-vowels', action='store', type=bool, default=True,
        help="Indicate whether you want to merge vowels in automatic segmentation.")
parser.add_argument('--expand-nasals', action='store', type=bool, default=False,
        help="Indicate whether you want to display nasals in vowels as an extra segment " +\
                'in automatic segmentation.')
parser.add_argument('-v', '--verbose', action='store_true', default=False, 
        help="Trigger verbose output.")

# subparsers
subparsers = parser.add_subparsers(dest="subcommand")

# currently, we define two programs for examples
parser_pairwise = subparsers.add_parser('pairwise')
parser_multiple = subparsers.add_parser('multiple')
parser_lexstat = subparsers.add_parser('lexstat')
parser_wordlist = subparsers.add_parser('wordlist')
parser_settings = subparsers.add_parser('settings')
parser_alignments = subparsers.add_parser('alignments')


for p in [parser_pairwise, parser_multiple, parser_lexstat, parser_alignments]:
    # we take care of pairwise first
    p.add_argument('-i', '--input-file', action='store', 
            help="Select your input file.")
    p.add_argument('-o', '--output-file', action='store', default='output-file',
            help="Select the name of the output file.")
    p.add_argument('--output', choices=["file", "stdout"], 
            default="stdout",  
            help="Select whether you want to have the output written to file or to screen.")
    p.add_argument('--factor', action='store', type=float, default=0.3,
            help="Select the factor for prosodic strings in SCA alignments.")
    p.add_argument('--gop', action='store', type=float, default=-1,
            help="Select your gap opening penalty.")
    p.add_argument('--scale', action='store', type=float, default=0.5,
            help="Select the gap extension scale.")
    p.add_argument('--restricted-chars', action='store', default='T_',
            help="Select the characters in the sound classes which serve to mark word " +\
                "or morpheme boundaries.")

# settings
parser_settings.add_argument('-p', '--params', action='store', nargs='+', 
        help="Indicate for which parameters you want to see the lingpy default settings.")

# we take care of pairwise first
parser_pairwise.add_argument('-s', '--strings', action='store', nargs=2,
        type=str,
        help="Type in the input sequences you want to align.")
parser_pairwise.add_argument('-m', '--mode',  action='store', default='global', 
        choices=['global', 'local', 'overlap', 'dialign'],
        help="Select the alignment mode.")
parser_pairwise.add_argument('-d', '--distance', action='store_true', default=False, 
        help="Choose whether you want distances or similarities to be reported.")
parser_pairwise.add_argument('--method', action='store', choices=['sca', 'basic'],
        default='basic', help="Select the basic method you want to use.")

# argparse for multiple
parser_multiple.add_argument('-s', '--strings', action='store', nargs="+",
        type=str,
        help="Type in the input sequences you want to align.")
parser_multiple.add_argument('-m', '--mode',  action='store', default='global', 
        choices=['global', 'overlap', 'dialign'], 
        help="Select the alignment mode.")
parser_multiple.add_argument('--method', action='store', choices=['sca', 'basic'],
        default='basic', help="Select the basic method you want to use.")
parser_multiple.add_argument('--tree-calc', action='store', choices=['ugpma', 'neighbor'],
        default='upgma', help="Select the tree cluster method you want to use for the guide tree.")
parser_multiple.add_argument('--gap-weight', action='store', type=float,
        default=0.5, help="Select the tree cluster method you want to use for the guide tree.")
parser_multiple.add_argument('--format', action='store', default='msa', 
        choices=['msa', 'html', 'tex', 'psa'],
        help="Select the output format.")
parser_multiple.add_argument('--align-method', action='store', default='library', 
        choices=['library', 'progressive'], 
        help="Select how you want to align, based on a library (T-COFFEE), or " +\
                "in a classical progressive way.")
parser_multiple.add_argument('--iteration', action='store', nargs="+", default=[],
        choices=['orphans', 'clusters', 'similar-gap-sites', 'all-sequences'],
        help="Select method for postprocessing.")
parser_multiple.add_argument('--swap-check', action='store_true', default=False, 
        help="Indicate whether you want to search for swapped sites automatically.")

# args for alms (alignment of cognates in wordlists)
parser_alignments.add_argument('-m', '--mode',  action='store', default='global', 
        choices=['global', 'overlap', 'dialign'], 
        help="Select the alignment mode.")
parser_alignments.add_argument('--method', action='store', choices=['sca', 'basic'],
        default='basic', help="Select the basic method you want to use.")
parser_alignments.add_argument('--tree-calc', action='store', choices=['ugpma', 'neighbor'],
        default='upgma', help="Select the tree cluster method you want to use for the guide tree.")
parser_alignments.add_argument('--gap-weight', action='store', type=float,
        default=0.5, help="Select the tree cluster method you want to use for the guide tree.")
parser_alignments.add_argument('--format', action='store', default='tsv', 
        choices=['msa', 'html', 'tsv'],
        help="Select the output format.")
parser_alignments.add_argument('--align-method', action='store', default='library', 
        choices=['library', 'progressive'], 
        help="Select how you want to align, based on a library (T-COFFEE), or " +\
                "in a classical progressive way.")
parser_alignments.add_argument('--iterate', action='store_true', default=False,
        help="Choose whether you want to make postprocessing using iteration methods.")
parser_alignments.add_argument('-c','--cognate-identifier', action='store',
        default='lingpyid', 
        help='Select the name for the column with the cognate judgments.')
parser_alignments.add_argument('--alignment-identifier', action='store',
        default='alignment', 
        help='Select the name for the column with the alignments.')
parser_alignments.add_argument('--use-logodds', action='store_true',
        default=False, 
        help='Use logodds from LexStat to calculate the alignments.')

# args for lexstat
parser_lexstat.add_argument('--method', action='store', default='sca',
        choices=['lexstat', 'sca','edit-dist','turchin'], 
        help='Select the method you want to use for cognate detection.')
parser_lexstat.add_argument('--runs', action='store', type=int, default=1000,
        help='Select the number of iterations for the permutation test.')
parser_lexstat.add_argument('--check', action='store_true', default=False, 
        help='Specify whether you want to carry out a check before running your analysis.')
parser_lexstat.add_argument('--scoring-threshold', action='store',
        default=lingpy.settings.rcParams['lexstat_scoring_threshold'],
        help='Select the threshold to be used when creating the scorer.')
parser_lexstat.add_argument('-t', '--threshold', action='store', type=float,
        default=0.45,
        help='Select the threshold to be used when creating the scorer.')
parser_lexstat.add_argument('-c','--cognate-identifier', action='store',
        default='lingpyid', 
        help='Select the name for the column with the cognate judgments.')
parser_lexstat.add_argument('--cluster-method', action='store',
        default='upgma', 
        choices=['upgma', 'single', 'complete', 'mcl', 'ward',
        'link_clustering'],
        help='Select the name of the cluster method you want to use.')
parser_lexstat.add_argument('--ratio', action='store', nargs=2, type=int,
        default=[2,1],
        help='Select the ratio between the logodds and the sound class scorer.')


def alignments(args):
    """
    Carry out alignment analysis of a wordlist file with readily detected cognates.

    """
    if args.input_file:
        alms = lingpy.align.sca.Alignments(args.input_file, 
                merge_vowels=args.merge_vowels,
                expand_nasals=args.expand_nasals,
                ref=args.cognate_identifier,
                alignment=args.alignment_identifier)
        if args.use_logodds:
            lex = lingpy.compare.lexstat.LexStat(args.input_file)
            if not hasattr(lex, 'cscorer'):
                raise ValueError("No scorer has been submitted along with your file.")
            scoredict = lex.cscorer
        else:
            scoredict = False
        alms.align(method=args.align_method, gop=args.gop, scale=args.scale,
                factor=args.factor, iteration=args.iterate, model=args.model,
                tree_calc=args.tree_calc,
                restricted_chars=args.restricted_chars, mode=args.mode)
        alms.output(args.format, filename=args.output_file, ignore='all')

        return alms.msa[args.cognate_identifier]

def settings(args):
    out = []
    for param in args.params:
        if param in lingpy.settings.rcParams:
            print('{0:20} : {1}'.format(param, repr(lingpy.settings.rc(param))))
        else:
            print('{0:20} : PARAMETER NOT FOUND'.format(param))
        out += [[param, lingpy.settings.rc(param) if param in
            lingpy.settings.rcParams else 'PARAMETER NOT FOUND']]
    return out

def lexstat(args):
    if args.input_file:
        lex = lingpy.compare.lexstat.LexStat(args.input_file, 
                merge_vowels=args.merge_vowels,
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
                    print("[i] Scorer has already been calculated. Skipping re-calculation.")
        lex.cluster(args.method, threshold=args.threshold,
                cluster_method=args.cluster_method, scale=args.scale,
                factor=args.factor, restricted_chars=args.restricted_chars,
                gop=args.gop, ref=args.cognate_identifier)
        if args.verbose: print("[i] Clustered the data.")
        lex.output('tsv', filename=args.output_file)
        return len(lex.get_etymdict(ref=args.cognate_identifier))

def pairwise(args):
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
    make_out = lambda x, y, z: '\t'.join(x)+'\n'+'\t'.join(y)+'\n{0:.2f}'.format(z)
    if args.strings:
        if args.method == 'basic':
            almA, almB, sim = lingpy.align.pairwise.pw_align(
                    args.strings[0], 
                    args.strings[1],
                    **args.__dict__
                    )
        else:
            pair = lingpy.align.pairwise.Pairwise(*args.strings, merge_vowels=args.merge_vowels,
                    expand_nasals=args.expand_nasals)
            pair.align(**args.__dict__)
            almA, almB, sim = pair.alignments[0]
        output = make_out(almA, almB, sim)
        if args.output == 'file':
            lingpy.util.write_text_file(args.output_file, output)
        else:
            print(output)
        return almA, almB, sim
    elif args.input_file:
        pairs = lingpy.align.sca.PSA(args.input_file)
        output = ''
        if args.method == 'basic':
            for i, (seqA, seqB) in enumerate(pairs.tokens):
                almA, almB, sim = lingpy.align.pairwise.pw_align(
                        ''.join(seqA),
                        ''.join(seqB)
                        ,
                        **args.__dict__)
                pairs.alignments[i] = [almA, almB, sim]
                output += make_out(almA, almB, sim)+'\n'
        elif args.method == 'sca':
            pairs.align(**args.__dict__)
            for almA, almB, sim in pairs.alignments:
                output += make_out(almA, almB, sim)+'\n'
        if args.output == 'file':
            pairs.output('psa', filename=args.output_file)
        else:
            print(output)
        return pairs.alignments

def multiple(args):
    """
    Multiple alignment console interface for LingPy.

    Notes
    -----

    """
    def align_msa(msa, args):
        if args.method == 'library':
            msa.lib_align(**args.__dict__)
        else:
            msa.prog_align(**args.__dict__)
        if 'orphans' in args.iteration: msa.iterate_orphans()
        if 'clusters' in args.iteration: msa.iterate_clusters()
        if 'similar-gap-sites' in args.iteration: msa.iterate_similar_gap_sites()
        if 'all-sequences' in args.iteration: msa.iterate_all_sequences()
        if args.swap_check: msa.swap_check()

    if args.strings:
        if args.method == 'basic':
            alms = lingpy.align.multiple.mult_align(
                    args.strings,
                    pprint=False,
                    gop=args.gop,
                    tree_calc=args.tree_calc,
                    scoredict=False,
                    scale=args.scale,
                    )
        elif args.method == 'sca':
            msa = lingpy.align.multiple.Multiple(args.strings, merge_vowels=args.merge_vowels,
                    expand_nasals=args.expand_nasals, model=args.model)
            align_msa(msa, args)
            alms = msa.alm_matrix
        
        output = '\n'.join(['\t'.join(seq) for seq in alms])
        if args.output == 'file':
            lingpy.util.write_text_file(args.output_file, output)
        else:
            print(output)
        return alms
    elif args.input_file:
        msa = lingpy.align.sca.MSA(args.input_file, merge_vowels=args.merge_vowels,
                    expand_nasals=args.expand_nasals, model=args.model)
        if args.method == 'basic':
            seqs = [str(''.join(t)) for t in msa.tokens]
            print(seqs)
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
        if args.output == 'file':
            msa.output(args.format, filename=args.output_file)
        else:
            for taxon, alm in zip(msa.taxa, msa.alm_matrix):
                print('{0}\t{1}'.format(taxon, '\t'.join(alm)))
        return msa.taxa, msa.alm_matrix
    
def main():
    """
    Basic command line interface function for LingPy

    Notes
    -----
    Major options for main in global are:

    * schema (ipa, asjp): select the alphabet underlying the input strings
    * model (sca, dolgo, asjp): selct the major sound-class model to be used

    
    """
    # parse arguments
    args = parser.parse_args()
    
    # make up for basic configs
    lingpy.rc(schema=args.schema)
    lingpy.rc(model=lingpy.settings.rcParams[args.model])
    
    if args.subcommand == "pairwise":
        pairwise(args)
    if args.subcommand == 'multiple':
        multiple(args)
    if args.subcommand == 'settings':
        settings(args)
    if args.subcommand == 'lexstat':
        lexstat(args)
    if args.subcommand == 'alignments':
        alignments(args)


