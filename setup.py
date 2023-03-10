from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='cghla',
      version='0.2.3',
      description='Tools for finding HLA alleles for WGS/RNAseq samples',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Marcus Breese',
      author_email='marcus@breese.com',
      url='http://github.com/compgen-io/cghla/',
      packages=['cghla'],
      scripts=['bin/cghla'],
      python_requires='>=3.9',

     )

