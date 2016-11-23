from setuptools import find_packages, setup

with open('README.rst', 'r') as f:
    long_description = f.read()

setup(
    name='nexthiv',
    version='0.0.1',
    author='Simon Frost',
    author_email='sdwfrost@gmail.com',
    packages=find_packages(exclude=['docs']),
    url='https://github.com/sdwfrost/nexthiv',
    license='MIT',
    description='Processing code for nextHIV',
    long_description=long_description,
    install_requires=[
    ]
)
