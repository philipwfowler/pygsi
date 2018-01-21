from setuptools import setup

setup(
    install_requires=[
        "numpy >= 1.13",
        "requests >= 2.8",
        "pandas >= 0.21.0",
        "tqdm >= 4.19"
    ],
    name='pygsi',
    version='0.1.0',
    url='https://github.com/philipwfowler/pygsi',
    author='Philip W Fowler',
    packages=['pygsi'],
    license='MIT',
    long_description=open('README.md').read(),
)
