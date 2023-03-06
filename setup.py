from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='cgmhc',
      version='0.1.1',
      description='Tools for finding MHC alleles for WGS/RNAseq samples',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Marcus Breese',
      author_email='marcus@breese.com',
      url='http://github.com/compgen-io/cgmhc/',
      packages=['cgmhc'],
      scripts=['bin/cgmhc'],
      python_requires='>=3.9',

     )

