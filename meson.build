project('rfpy', 'c',
  version : '0.2.0',
  license: 'MIT',
  meson_version: '>=0.64.0',
)

py_mod = import('python')
py = py_mod.find_installation()
py_dep = py.dependency()

py.install_sources(
  'rfpy/rfdata.py',
  'rfpy/plotting.py',
  'rfpy/hk.py',
  'rfpy/harmonics.py',
  'rfpy/binning.py',
  'rfpy/utils.py',
  'rfpy/__init__.py',
  'rfpy/examples/data/*',
  'rfpy/scripts/*',
  pure: false,
  subdir: 'rfpy'
  )
