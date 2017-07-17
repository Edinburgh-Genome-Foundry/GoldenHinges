import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('goldenhinges/version.py').read()) # loads __version__

setup(
	name='goldenhinges',
    version=__version__,
    author='Zulko',
    description='',
    long_description=open('README.rst').read(),
    license='see LICENSE.txt',
    keywords="golden gate assembly overhangs constraint programming",
    packages= find_packages(exclude='docs'),
    install_requires=['numpy', 'Numberjack', 'Biopython']
)
