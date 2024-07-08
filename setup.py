import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Pike",
    version="0.1.0",
    author="D.V. Krivonos",
    author_email="danil01060106@gmail.com",
    description="A tool for Oxford Nanopore amplicon metagenomics",
    long_description="Pike",
    long_description_content_type="",
    url="https://github.com/DanilKrivonos/Pike",
    project_urls={
        "Bug Tracker": "https://github.com/DanilKrivonos/Pike",
    },
    package_data={

    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    include_package_data=True,
    packages=['pike', 'Pike.src', 'Pike.example_data', 'get_taxonomy'],
    install_requires=[
        'numpy',
        'scikit-bio>=0.5.9'
        'biopython>=1.83',
        'pandas>=2.2.0',
        'scipy>=1.10.1',
        'scikit-learn>=1.4.1',
        'hdbscan>=0.8.33',
        'umap-learn>=0.5.5',
        ''
    ],
    entry_points={
        'console_scripts': [
            'pike=Pike.pike:main',
            'get_taxonomy=Pike.get_taxonomy:main'
        ]
    }
)