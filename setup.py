import setuptools
from distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='snpy',
    version='0.1.0',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['snpy',],
    package_dir={'snpy': 'snpy'},
    python_requires='>=3.6',
)
