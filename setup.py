import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('goldenhinges/version.py').read()) # loads __version__

setup(
	name='goldenhinges',
    version=__version__,
    author='Zulko',
    description='DNA overhangs design for Golden Gate etc.',
    url='https://github.com/Edinburgh-Genome-Foundry/GoldenHinges',
    long_description=open('pypi-readme.rst').read(),
    license='see LICENSE.txt',
    keywords="DNA assembly overhangs constraint-programming synthetic-biology",
    packages= find_packages(exclude='docs'),
    install_requires=['numpy', 'Numberjack', 'Biopython', 'networkx',
                      'proglog', 'dna_features_viewer', 'tqdm', 'dnachisel',
                      'geneblocks']
)
