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
    zip_safe=False,
    license='MIT',
    description='Processing code for nextHIV',
    long_description=long_description,
    install_requires=[
        'rethinkdb>=2.3',
        'PyYAML>=3.12',
        'biopython>=1.68',
        'pandas>=0.19.1',
        'rpy2>=2.8.4'
    ],
    dependency_links=[
        'git+ssh://git@github.com/veg/bioext.git@master'
    ],
    scripts=[
        'bin/nexthiv'
    ],
    data_files=[
        ('nexthiv/data',[
            'nexthiv/data/Scores_PI.txt',
            'nexthiv/data/Scores_NRTI.txt',
            'nexthiv/data/Scores_NNRTI.txt',
            'nexthiv/data/combinationScores_PI.txt',
            'nexthiv/data/combinationScores_NRTI.txt',
            'nexthiv/data/combinationScores_NNRTI.txt',
            'nexthiv/data/hiv_refs_prrt_trim.bam',
            'nexthiv/data/hiv_refs_prrt_trim.bam.bai',
            'nexthiv/data/SDRM_2009.txt'
        ])
    ]
)
