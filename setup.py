from setuptools import setup, find_packages

setup(
    name='g3py',
    version='0.1',
    description='A Python Library for Gamma Ray Analysis',
    author='mohankarthikg',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'astropy'  
    ],
)
