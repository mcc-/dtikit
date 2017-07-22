import os
import logging
import shelve

_config_dir = os.path.expanduser('~'+os.sep+'.dtikit'+os.sep)
if not os.path.isdir(__config_dir):
    os.mkdir(__config_dir)
_config = shelve.open(filename = __config_dir+os.sep+'config.yaml', protocol=2)
if 'last_folder' not in _config:
    _config['last_folder'] = None
if 'log_file' not in _config:
    _config['log_file'] = _config_dir+'dtikit.log'
logging.basicConfig(filename=_config['log_file'])
