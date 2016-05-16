# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
import sys
from contextlib import contextmanager

from six import StringIO

from lingpy.cli import main
from lingpy.tests.util import test_data, WithTempDir


@contextmanager
def capture(*args):
    out, sys.stdout = sys.stdout, StringIO()
    main(*args)
    sys.stdout.seek(0)
    yield sys.stdout.read()
    sys.stdout = out


class Tests(WithTempDir):
    def run_cli(self, *args):
        if len(args) == 1:
            args = args[0].split()
        with capture(*args) as output:
            return output

    def test_wordlist(self):
        output = self.run_cli(
            'wordlist -i {0} --stats --calculate diversity'.format(test_data('KSL.qlc')))
        self.assertIn('Height:  200', output)

    def test_alignments(self):
        tmp = self.tmp_path('alignments')

        def cmd(i, rem=''):
            return 'alignments -i {0} -c cogid -o {1} {2}'.format(i, tmp.as_posix(), rem)

        self.run_cli(cmd(test_data('KSL.qlc')))
        self.run_cli(cmd(test_data('KSL3.qlc'), ' --format html --use-logodds'))

    def test_multiple(self):
        # first test, align string, no output, no input
        output = self.run_cli('multiple -s woldemort waldemar walter')
        self.assertIn('w\ta\tl\tt\te\t-\t-\tr\t-', output)

        # second test, test output as file, no input, vary method as sca
        mlt = main('multiple', '-s', 'woldemort', 'waldemar',
                                  'walter', '--method', 'sca', '--output-file',
                                  self.tmp_path('out.msa').as_posix())
        # third test, test output and input
        # second test, test output as file, no input, vary method as sca
        mlt = main('multiple', '-i', test_data('harryp.msa'),
                                  '--method', 'sca', '--output-file',
                                  self.tmp_path('out2.msa').as_posix(), '--align-method',
                                  'library')
        assert len(mlt[0]) == 7

        # fourth test, test output and input with method=basic
        mlt = main('multiple', '-i', test_data('harryp.msa'),
                                  '--method', 'basic', '--output-file',
                                  self.tmp_path('out2.msa').as_posix())
        assert len(mlt[0]) == 7
        assert len([x for x in mlt[1][-1] if x != '-']) == 4

    def test_pairwise(self):
        # first test, align string, no output, no input
        output = self.run_cli('pairwise -s woldemort waldemar')
        self.assertEqual(
            [line.split('\t') for line in output.split('\n')][:2],
            [
                ['w', 'o', 'l', 'd', 'e', 'm', 'o', 'r', 't'],
                ['w', 'a', 'l', 'd', 'e', 'm', 'a', 'r', '-'],
            ])

        # second test, test output as file, no input, vary method as sca
        tmp = self.tmp_path('test1')
        self.run_cli(
            'pairwise -s woldemort waldemar --method sca -o {0} --distance'.format(
                tmp.as_posix()))
        assert tmp.exists()

        # third test, test output and input
        # second test, test output as file, no input, vary method as sca
        tmp = self.tmp_path('test2')
        self.run_cli('pairwise -i {0} --method sca -o {1} -m overlap'.format(
            test_data('harry_potter.psa'), tmp.as_posix()))
        #
        # FIXME: It should not be the case that an output file is soecified, but the
        # output is actually written to a different file!
        #
        assert tmp.parent.joinpath(tmp.name + '.psa').exists()

        # fourth test, test output and input with method=basic
        tmp = self.tmp_path('test3')
        self.run_cli('pairwise -i {0} --method basic -o {1}'.format(
            test_data('harry_potter.psa'), tmp.as_posix()))
        assert tmp.parent.joinpath(tmp.name + '.psa').exists()

    def test_settings(self):
        output = self.run_cli('settings -p lexstat_threshold lexstat_runs')
        self.assertNotIn('gop', output)
        self.assertIn('lexstat_runs', output)
        self.assertIn('lexstat_threshold', output)

    def test_lexstat(self):
        cogs = main('lexstat', '-i',
                                  test_data('KSL.qlc'),
                                  '--output-file', self.tmp_path('lexstat').as_posix())
        assert cogs == 1080
