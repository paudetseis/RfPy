project('rfpy', 'c',
  version : '0.2.1',
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
  'rfpy/ccp.py',
  'rfpy/__init__.py',
  'rfpy/scripts/rfpy_calc.py',
  'rfpy/scripts/rfpy_recalc.py',
  'rfpy/scripts/rfpy_plot.py',
  'rfpy/scripts/rfpy_hk.py',
  'rfpy/scripts/rfpy_harmonics.py',
  'rfpy/scripts/rfpy_ccp.py',
  'rfpy/scripts/__init__.py',
  pure: false,
  preserve_path: true,
  )
  