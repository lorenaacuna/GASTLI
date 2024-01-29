

project('gastli', 'c',
  version : '0.1',
  license: 'BSD-3',
  meson_version: '>=0.64.0',
    default_options: [
    'buildtype=debugoptimized',
    # TODO: the below -Wno flags are all needed to silence warnings in
    # f2py-generated code. This should be fixed in f2py itself.
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types',
    'fortran_args=-Wno-conversion',  # silence "conversion from REAL(8) to INTEGER(4)"
    'fortran_std=legacy',
    'warning_level=3']
)

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

#https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()
message(py3.full_path())
message(py3.get_install_dir())

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
  ['-c', 'import os; os.chdir(".."); from numpy import f2py; print(f2py.get_include())'],
  check : true
).stdout().strip()

inc_dirs = include_directories(incdir_numpy, incdir_f2py)

# Platform detection to set more flags for Windows systems
is_windows = host_machine.system() == 'windows'
is_mac = host_machine.system() == 'darwin'
is_mingw = is_windows and cc.get_id() == 'gcc'

message(host_machine.system())

# Check the python headers
# cc.check_header('Python.h', dependencies: [py3_dep], required: true)

# Unlike distutils, meson doesn't yet include some of the f2py stuff
fortranobject_c = incdir_f2py / 'fortranobject.c'

fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  dependencies: py3_dep,
  include_directories: [incdir_numpy, incdir_f2py])

fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [incdir_numpy, incdir_f2py])

# We create a custom target to compile and make the module
dimensions_source = custom_target('dimensionsmodule.c',
  input : ['src/mod_dimensions.f95'],
  output : ['dimensionsmodule.c', 'dimensions-f2pywrappers2.f90'],
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'dimensions', '--lower']
  )

py3.extension_module('dimensions',
  ['src/mod_dimensions.f95', dimensions_source],
  include_directories: inc_dirs,
  dependencies : [py3_dep, fortranobject_dep],
  subdir: 'gastli',   # if you want to install this part as a submodule
  install : true
)

py3.install_sources([],
  pure: false,
  subdir: 'gastli'
)

