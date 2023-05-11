import inspect
import os
import sys
from pathlib import Path
import cmake_build_extension
import setuptools
from distutils.util import convert_path

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
long_description = (this_directory / "README.md").read_text()

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
    description='python bead tracker',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Mriv31/PyTrackBead',
    license='BSD 3-Clause License',
    author='Martin Rieu',
    author_email='martin.rieu@physics.ox.ac.uk',
    cmdclass=dict(
        build_ext=cmake_build_extension.BuildExtension,
    ),
)
