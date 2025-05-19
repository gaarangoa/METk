from setuptools import setup, find_packages

setup(
    name='metk',
    version='1.0.1',
    author='Gustavo Arango',
    author_email='',
    description='Mutation Enrichment Toolkit',
    packages=find_packages(),
    install_requires=[
        # List your module's dependencies here
    ],
    entry_points={
        'console_scripts': [],
    },
    include_package_data=True,
    package_data={
        'metk': ['embed_doc']
    },
)