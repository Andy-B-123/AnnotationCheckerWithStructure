# setup.py

from setuptools import setup, find_packages

setup(
    name="AnnotationCheckerWithStructure",
    version="0.1.5",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        "biopython==1.79",
        "BCBio-GFF==0.7.1",
        "pandas==1.3.5",
        "matplotlib==3.5.1",
        "seaborn==0.11.2"
    ],
    entry_points={
        'console_scripts': [
            'AnnotationCheckerWithStructure=main:main',
        ],
    },
    author="Andreas Bachler",
    author_email="Andy.Bachler@example.com",
    description="A simple bioinformatics script.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/Andy-B-123/AnnotationCheckerWithStructure",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.11',
)
