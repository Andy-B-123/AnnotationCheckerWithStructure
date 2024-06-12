# setup.py

from setuptools import setup, find_packages

setup(
    name="AnnotationCheckerWithStructure",
    version="0.1.2",
    packages=find_packages(),
    install_requires=[
        "pandas>=2.0.3"
    ],
    entry_points={
        'console_scripts': [
            'AnnotationCheckerWithStructure=AnnotationCheckerWithStructure.main:main',
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
    python_requires='>=3.6',
)
