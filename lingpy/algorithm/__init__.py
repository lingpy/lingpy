try:
    from .calignx import *
except:
    print("[i] Loading C module failed, using pure Python module...")
    from .c_alignx import *
from .cluster import *
from .misc import *
