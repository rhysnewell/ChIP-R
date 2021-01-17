from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='ChIP-R',
    version='1.2.0',
    author='Rhys Newell, Mikael Boden, Alex Essebier',
    author_email='r.newell@uq.edu.au',
    packages=find_packages(),
    license='GPL-3.0',
    long_description_content_type="text/markdown",
    long_description=long_description,
    scripts=['bin/chipr'],
    install_requires = ['scipy', 'numpy'],
    entry_points={
        'console_scripts': [
        'chipr=chipr.__main__:main',
        'ChIP-R=chipr.__main__:main'
        ]
    },
    python_requires='>=3',
    description=('ChIP-R is a method for assessing the reproducibility of ' +
                 'replicated ChIP-seq type experiments. It incorporates the ' +
                 'rank product method, a novel thresholding methods, and the ' +
                 'use of peak fragmentation return the most reproducible peaks.'),
    keywords='ChIP-R',
    url='https://github.com/rhysnewell/ChIP-R',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
)

