import os

_PKG_DIR = os.path.expanduser(os.sep.join(['~','.dtikit']))
if not os.path.isdir(_PKG_DIR):
  os.mkdir(_PKG_DIR)
