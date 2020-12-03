import pytest

from lingpy.config import Config

CFG_NAME = '.config'


@pytest.fixture
def config_factory(tmppath):
    def make(**kw):
        return Config('test', config_dir=tmppath / CFG_NAME, **kw)
    return make


def test_new_config(config_factory, tmppath):
    assert not tmppath.joinpath(CFG_NAME).exists()
    cfg = config_factory()
    assert tmppath.joinpath(CFG_NAME).exists()


def test_existing_config(config_factory, tmppath):
    tmppath.joinpath(CFG_NAME).mkdir()
    tmppath.joinpath(CFG_NAME, 'test.ini').write_text(
        """\
[section]
option = 12
""", encoding='utf8')
    cfg = config_factory()
    assert cfg.get('section', 'option') == '12'


def test_default(config_factory):
    cfg = config_factory(default="""\
[section2]
option2 = 7
""")
    assert cfg.get('section2', 'option2') == '7'
