import ez_setup

ez_setup.use_setuptools()

from setuptools import setup, find_packages

version = {}
with open("goldenhinges/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="goldenhinges",
    version=version["__version__"],
    author="Zulko",
    description="DNA overhang design for Golden Gate etc.",
    url="https://github.com/Edinburgh-Genome-Foundry/GoldenHinges",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    keywords="DNA assembly overhangs constraint-programming synthetic-biology",
    packages=find_packages(exclude="docs"),
    install_requires=[
        "numpy",
        "Numberjack",
        "Biopython",
        "networkx",
        "proglog",
        "dna_features_viewer",
        "tqdm",
        "dnachisel",
        "geneblocks",
    ],
)
