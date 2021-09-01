from setuptools import find_packages
import setuptools

VERSION = '1.0.1'
DESCRIPTION = 'EpIntervene Simulation Package'
LONG_DESCRIPTION = 'EpIntervene package for simulating SIR and SEIR Epidemics with or without interventions'

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="epintervene",
    version=VERSION,
    author="Andrea Allen",
    author_email="andrea2allen@gmail.com",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=[],  # add any additional packages that
    # needs to be installed along with your package.

    keywords=['python', 'epintervene'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires=">=3.6",
)