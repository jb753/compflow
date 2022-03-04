import setuptools

from numpy.distutils.core import setup
from numpy.distutils.core import Extension

ext = [Extension(name = 'compflow_fort_from_Ma', sources = ['fortran/from_Ma.f']),
       Extension(name = 'compflow_fort_der_from_Ma', sources = ['fortran/der_from_Ma.f']),
       Extension(name = 'compflow_fort_to_Ma', sources = ['fortran/to_Ma.f'])]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="compflow",
    version="0.4.0",
    author="James Brind",
    author_email="jb753@cam.ac.uk",
    description="A library of fast compressible flow relations.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://jamesbrind.uk/compflow-docs/index.html",
    packages=setuptools.find_packages(),
    install_requires=[
        "numpy",
        "scipy"
        ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    keywords=['compressible','flow','aerodynamics','fluid','engineering'],
    ext_modules=ext,
)
