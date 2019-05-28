from setuptools import setup, find_packages

setup(
    name='chipr',
    version='1.0',
    author='Rhys Newell',
    suthor_email='r.newell@uq.net.au',
    packages=find_packages(),
    license='GPL-3.0',
    long_description=open('README.md').read(),
    scripts=['bin/chipr'],
    install_requires = ['scipy', 'numpy'],
    entry_points={
        'console_scripts': [
        'chipr=chipr.__main__:main'
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

