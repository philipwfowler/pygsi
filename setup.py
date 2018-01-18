from distutils.core import setup

setup(
    install_requires=[
        "numpy",
        "json",
        "requests",
        "pandas"
    ],
    name='pygsi',
    version='0.1.0',
    url='https://github.com/philipwfowler/pygsi',
    author='Philip W Fowler',
    packages=['pygsi'],
    license='MIT',
    long_description=open('README.md').read(),
)
