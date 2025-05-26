from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, Extension
import pybind11
import os

# Define the extension module
ext_modules = [
    Pybind11Extension(
        "spz._core",
        [
            "src/cc/load-spz.cc",
            "src/cc/splat-c-types.cc", 
            "src/cc/splat-types.cc",
            "python/spz_bindings.cpp"
        ],
        include_dirs=[
            pybind11.get_cmake_dir() + "/../../../include",
            "src/cc"
        ],
        libraries=["z"],  # zlib
        cxx_std=17,
    ),
]

setup(
    name="spz",
    version="1.1.0",
    author="Niantic",
    description="Python bindings for SPZ compressed 3D gaussian splats",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
    packages=["spz"],
    package_dir={"spz": "python/spz"},
    install_requires=[
        "numpy>=1.19.0",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
    ],
) 