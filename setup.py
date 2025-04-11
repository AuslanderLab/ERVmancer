from setuptools import setup, find_packages

setup(
    name="ervmancer",
    version="0.1.0",
    description="ERVmancer is a bioinformatics tool to quantify Human Endogenous Retrovirus (HERV) short read RNA sequencing expression data by aligning short reads to a curated subset of HERVs and then resolving ambiguity in alignment using a pre-computed HERV phylogenetic tree.",
    author="Andrew Patterson",
    author_email="Andrew.Patterson@pennmedicine.upenn.edu",
    url="https://github.com/AuslanderLab/ervmancer",
    packages=find_packages("src"),
    package_dir={"": "src"},
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "pandas",
        "tqdm",
        "regex",
        "setuptools"
    ],
    entry_points={
        'console_scripts': [
            'ervmancer=ervmancer.main:main',
        ],
    },
)
