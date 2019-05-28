from setuptools import setup, find_packages

setup(
    name='chipr',
    version='1.0',
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
)

