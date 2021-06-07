"""
Setup-Script for LingPy
"""
from setuptools import setup, find_packages


setup(
    name="lingpy",
    description="Python library for quantitative tasks in historical linguistics",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    author="Johann-Mattis List and Simon Greenhill and Robert Forkel",
    author_email="info@lingpy.org",
    version="2.6.8",
    packages=find_packages(
        where="src",
        exclude=[
            "lingpy._plugins",
            "_plugins",
            "*._plugins",
            "_plugins.*",
            "*._plugins.*",
            "build",
            "private",
            "lingpy.egg-info",
            "dist",
            "lib",
        ],
    ),
    package_dir={"": "src"},
    install_requires=[
        "numpy",
        "appdirs",
        "networkx>=2.3",
        "tqdm",
        'csvw>=1.5.6"',
        'clldutils>=2.8.0',
        'pycldf>=1.7.0',
    ],
    extras_require={
        "borrowing": ["matplotlib", "scipy"],
        "cluster": ["python-igraph", "scikit-learn"],
        "test": ["pytest", "coverage", "pytest-mock", "pytest-cov"],
        "dev": ["wheel", "twine", "sphinx", "tox"],
    },
    entry_points={"console_scripts": ["lingpy=lingpy.cli:main"]},
    keywords=[
        "historical linguistics",
        "sequence alignment",
        "computational linguistics",
        "dialectology",
        "cognate detection",
    ],
    classifiers=[
        "Natural Language :: English",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering",
        "Topic :: Text Processing :: Linguistic",
    ],
    url="http://lingpy.org",
    license="gpl-3.0",
    platforms=["unix", "linux", "windows"],
    include_package_data=True,
    exclude_package_data={},
)
