import sys
from cx_Freeze import setup, Executable

base = None
if sys.platform == 'win32':
    base = 'Win32GUI'

options = {
    'build_exe': {
        'includes': 'atexit'
    }
}

executables = [
    Executable('graphlasso-covarience.py', base=base)
]

setup(name='gene_graph',
      version='0.1',
      description='Takes single cell RNA-sequencing and builds graph networks.',
      options=options,
      executables=executables
      )
