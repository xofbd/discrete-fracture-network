from setuptools import setup

setup(
    name='dfn',
    version='1.0.1',
    packages=['dfn'],
    author='Don Bruce Fox',
    author_email='dfox09@gmail.com',
    description='an analytical thermohydraulic model for discretely fractured\
    geothermal reservoirs',
    install_requires=['networkx', 'numpy', 'scipy'],
    test_suite='nose.collector',
    test_require=['nose'],
    license='MIT'
)
