from setuptools import setup, find_packages

setup(
    name="ervromancer",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        'tqdm',
    ],
    entry_points={
        'console_scripts': [
            'ervromancer=ervromancer.main:main',
        ],
    },
    python_requires=">=3.6",
)
