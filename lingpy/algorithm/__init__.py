try:
    from .alignx import *
except:
    print("[i] Loading C module failed, using pure Python module...")
    from ._alignx import *
from .cluster import *
from .misc import *
