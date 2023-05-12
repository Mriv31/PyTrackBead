import inspect
import os
import sys
from pathlib import Path
import cmake_build_extension
import setuptools
from distutils.util import convert_path, dir_util

from wheel.pep425tags import get_supported

if dir_util.is_file('dist/PyTrackBead-*.whl'):
    print('Found built wheel, installing...')
    whl_path = dir_util.get_path('dist/PyTrackBead-*.whl')
    whl_info = wheel.Wheel(whl_path)
    if whl_info.supported(version.LooseVersion(sys.version), sys.platform):
        pip install whl_path
    else:
        print('Wheel is not compatible with this Python version or architecture')

init_py = inspect.cleandoc(
    """
    import cmake_build_extension

    with cmake_build_extension.build_extension_env():
        from .bindings import *

    from .version import __version__

    """
)

CIBW_CMAKE_OPTIONS = [] # type: ignore
if "CIBUILDWHEEL" in os.environ and os.environ["CIBUILDWHEEL"] == "1":

    if sys.platform == "linux":
        CIBW_CMAKE_OPTIONS += ["-DCMAKE_INSTALL_LIBDIR=lib"]

version_path = convert_path('PyTrackBead/version.py')
version_dict = {}
with open(version_path) as version_file:
    exec(version_file.read(), version_dict)

from pathlib import Path
this_directory = Path(__file__).parent

setuptools.setup(
    ext_modules=[
        cmake_build_extension.CMakeExtension(
            name="Pybind11Bindings",
            install_prefix="PyTrackBead",
            cmake_depends_on=["pybind11"],
            write_top_level_init=init_py,
            source_dir=str(Path(__file__).parent.absolute()),
            cmake_configure_options=[
                f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
                "-DCALL_FROM_SETUP_PY:BOOL=ON",
                "-DBUILD_SHARED_LIBS:BOOL=OFF",
            ]
            + CIBW_CMAKE_OPTIONS,
        ),
    ],
    version=version_dict['__version__'],
    cmdclass=dict(
        build_ext=cmake_build_extension.BuildExtension,
    ),
)
