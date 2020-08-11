import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="compflow", # Replace with your own username
    version="0.1.0",
    author="James Brind",
    author_email="brind.james@gmail.com",
    description="Compressible flow tables",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jb753/compflow",
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
)
