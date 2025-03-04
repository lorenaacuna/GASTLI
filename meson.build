

project('gastli', 'c',
  version : '1.0.0',
  license: 'BSD-3',
  meson_version: '>=0.64.0',
    default_options: [
    'buildtype=debugoptimized',
    #'buildtype=release',
    # TODO: the below -Wno flags are all needed to silence warnings in
    # f2py-generated code. This should be fixed in f2py itself.
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types',
    'fortran_args=-Wno-conversion',   # silence "conversion from REAL(8) to INTEGER(4)"
    'fortran_std=legacy',
    'warning_level=3']
)


# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

# Optimisation flag for fortran
if fc.has_argument('-Ofast')
  add_project_arguments(['-Ofast', '-fno-finite-math-only'], language: 'fortran')
else
  add_project_arguments(['-O3'], language: 'fortran')
endif

# OPENMP flags for fortran
if fc.get_id() == 'nvidia_hpc'
  add_global_arguments('-mp', language : 'fortran')
  add_global_link_arguments('-mp', language : 'fortran')
elif fc.get_id() == 'nagfor'
  add_global_arguments('-openmp', language : 'fortran')
  add_global_link_arguments('-openmp', language : 'fortran')
elif fc.get_id() == 'gcc'
  add_global_arguments('-fopenmp', language : 'fortran')
  add_global_link_arguments('-fopenmp', language : 'fortran')
elif fc.get_id() == 'intel' or fc.get_id() == 'intel-llvm'
  add_global_arguments('-qopenmp', language : 'fortran')
  add_global_link_arguments('-qopenmp', language : 'fortran')
elif fc.get_id() == 'intel-cl' or fc.get_id() == 'intel-llvm-cl'
  add_global_arguments('/Qopenmp', language : 'fortran')
endif


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

## dimensions ##
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

## constants ##
constants_source = custom_target('constantsmodule.c',
  input : ['src/mod_constants.f95'],
  output : ['constantsmodule.c', 'constants-f2pywrappers2.f90'],
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'constants', '--lower']
  )

py3.extension_module('constants',
  ['src/mod_constants.f95', constants_source],
  include_directories: inc_dirs,
  dependencies : [py3_dep, fortranobject_dep],
  subdir: 'gastli',   # if you want to install this part as a submodule
  install : true
)


## fort_input ##
fort_input_source = custom_target('fortinputmodule.c',
  input : ['src/fort_input.f95'],
  output : ['fortinputmodule.c', 'fortinput-f2pywrappers2.f90'],
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'fortinput', '--lower']
  )


py3.extension_module('fortinput',
  ['src/fort_input.f95', fort_input_source],
  include_directories: inc_dirs,
  dependencies : [py3_dep, fortranobject_dep],
  subdir: 'gastli',   # if you want to install this part as a submodule
  install : true
)


## fort_spec ##
fort_spec_source = custom_target('fortspecmodule.c',
  input : ['src/fort_spec_func.f90','src/fort_spec.f90'],
  output : ['fortspecmodule.c', 'fortspec-f2pywrappers2.f90'],
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'fortspec', '--lower']
  )

py3.extension_module('fortspec',
  ['src/fort_spec_func.f90','src/fort_spec.f90', fort_spec_source],
  include_directories: inc_dirs,
  dependencies : [py3_dep, fortranobject_dep],
  subdir: 'gastli',   # if you want to install this part as a submodule
  install : true
)


## GASTLI_interior ##
gastli_interior_source = custom_target('gastli_interiormodule.c',
  input : ['src/mod_funcmazevet.f95',
           'src/mod_mazevet.f95',
           'src/mod_constants.f95',
           'src/mod_dimensions.f95',
           'src/mod_parameters.f95',
           'src/mod_functions.f95',
           'src/mod_modEOS.f95',
           'src/mod_subroutines.f95',
           'src/GASTLI_interior.f95'],
  output : ['gastli_interiormodule.c', 'gastli_interior-f2pywrappers2.f90'],
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'gastli_interior', '--lower']
  )

if cc.get_id() == 'clang'
  # Meson's OpenMP dependency does not work on mac yet
  omp = []
  link_args_parallel = ['-lgomp']
else
  omp = dependency('openmp')
  link_args_parallel = ['-fopenmp','-lgomp']
endif



py3.extension_module('gastli_interior',
  ['src/mod_funcmazevet.f95',
   'src/mod_mazevet.f95',
   'src/mod_constants.f95',
   'src/mod_dimensions.f95',
   'src/mod_parameters.f95',
   'src/mod_functions.f95',
   'src/mod_modEOS.f95',
   'src/mod_subroutines.f95', 'src/GASTLI_interior.f95', gastli_interior_source],
  include_directories: inc_dirs,
  dependencies : [py3_dep, fortranobject_dep, omp],
  link_args : link_args_parallel,
  subdir: 'gastli',   # if you want to install this part as a submodule
  install : true
)



py3.install_sources(['src/GASTLI.py',
                     'src/Guillot10.py',
                     'src/radius_hse.py',
                     'src/atm_models_interp.py',
                     'src/Coupling.py',
                     'src/Thermal_evolution.py',
                     'src/water_curves.py'],
  pure: false,
  subdir: 'gastli'
)


