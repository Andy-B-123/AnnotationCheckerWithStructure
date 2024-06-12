# setup.py

from setuptools import setup, find_packages

setup(
    name="bioinformatics_script",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas>=2.0.3"
    ],
    entry_points={
        'console_scripts': [
            'bioinformatics_script=bioinformatics_script.main:main',
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
