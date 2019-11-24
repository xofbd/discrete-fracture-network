from setuptools import setup

from dfn import __version__

with open('requirements.txt', 'r') as f:
    requirements = [line.strip() for line in f]


setup(
    name='dfn',
    version=__version__,
    packages=['dfn'],
    author='Don Bruce Fox',
    author_email='dfox09@gmail.com',
    description='an analytical thermohydraulic model for discretely fractured\
    geothermal reservoirs',
    install_requires=requirements,
    test_suite='nose.collector',
    test_require=['nose'],
    license='MIT'
)
