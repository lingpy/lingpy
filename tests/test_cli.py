import shlex

from lingpy.cli import main


def run(capsys, *args):
    if len(args) == 1:
        args = shlex.split(args[0])
    main(*args)
    return capsys.readouterr().out


def test_wordlist(capsys, test_data):
    output = run(
        capsys, 'wordlist -i {0} --stats --calculate diversity'.format(test_data / 'KSL.qlc'))
    assert 'Height:  200' in output


def test_alignments(tmp_path, capsys, test_data):
    def cmd(i, rem=''):
        return 'alignments -i {0} -c cogid -o {1} {2}'.format(i, tmp_path / 'alignments', rem)

    run(capsys, cmd(test_data / 'KSL.qlc'))
    run(capsys, cmd(test_data / 'KSL3.qlc', ' --format html --use-logodds'))


def test_ortho_profile(tmp_path, capsys, test_data):
    run(capsys,
        'profile -i ' + str(test_data / 'KSL.qlc') + ' --column ipa -o ' + str(tmp_path / 'ortho'))


def test_multiple(capsys, test_data, tmp_path):
    # first test, align string, no output, no input
    output = run(capsys, 'multiple -s woldemort waldemar walter')
    assert 'w\ta\tl\tt\te\t-\t-\tr\t-' in output

    # second test, test output as file, no input, vary method as sca
    _ = run(capsys, 'multiple', '-s', 'woldemort', 'waldemar',
             'walter', '--method', 'sca', '--output-file', str(tmp_path / 'out.msa'))

    # third test, test output and input
    # second test, test output as file, no input, vary method as sca
    mlt = main('multiple', '-i', str(test_data / 'harryp.msa'),
               '--method', 'sca', '--output-file',
               str(tmp_path / 'out2.msa'), '--align-method',
               'library')
    assert len(mlt[0]) == 7

    # fourth test, test output and input with method=basic
    mlt = main('multiple', '-i', str(test_data / 'harryp.msa'),
               '--method', 'basic', '--output-file',
               str(tmp_path / 'out2.msa'))
    assert len(mlt[0]) == 7
    assert len([x for x in mlt[1][-1] if x != '-']) == 4


def test_pairwise(capsys, tmp_path, test_data):
    # first test, align string, no output, no input
    output = run(capsys, 'pairwise -s woldemort waldemar')
    assert [line.split('\t') for line in output.split('\n')][:2] == \
        [
            ['w', 'o', 'l', 'd', 'e', 'm', 'o', 'r', 't'],
            ['w', 'a', 'l', 'd', 'e', 'm', 'a', 'r', '-'],
        ]

    # second test, test output as file, no input, vary method as sca
    tmp = tmp_path / 'test1'
    run(capsys, 'pairwise -s woldemort waldemar --method sca -o {0}'
                 ' --distance'.format(tmp))
    assert tmp.exists()

    # third test, test output and input
    # second test, test output as file, no input, vary method as sca
    tmp = tmp_path / 'test2'
    run(capsys, 'pairwise -i {0} --method sca -o {1} -m overlap'.format(
        test_data / 'harry_potter.psa', tmp))
    #
    # FIXME: It should not be the case that an output file is soecified, but the
    # output is actually written to a different file!
    #
    assert tmp.parent.joinpath(tmp.name + '.psa').exists()

    # fourth test, test output and input with method=basic
    tmp = tmp_path / 'test3'
    run(capsys, 'pairwise -i {0} --method basic -o {1}'.format(
        test_data / 'harry_potter.psa', tmp))
    assert tmp.parent.joinpath(tmp.name + '.psa').exists()


def test_settings(capsys):
    output = run(capsys, 'settings -p lexstat_threshold lexstat_runs')
    assert 'gop' not in output
    assert 'lexstat_runs' in output
    assert 'lexstat_threshold' in output


def test_lexstat(test_data, tmp_path):
    # TODO check with other python versions
    cogs = main('lexstat', '-i', str(test_data / 'KSL.qlc'),
                '--output-file', str(tmp_path / 'lexstat'))
    print(cogs)
    assert cogs in [1080, 1081]


def test_profile(test_data):
    prf = main('profile', '-i', str(test_data / 'KSL.qlc'), '--column', 'ipa')
    assert prf == 105
    prf = main('profile', '-i', str(test_data / 'KSL.qlc'), '--column', 'ipa',
               '--language', 'German', '--count')
    assert prf == 37
