from setuptools import setup, find_packages

setup(
    name="hvreglocator",
    version="0.1",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "hvreglocator=hvreglocator.main:main",
        ],
    },
    install_requires=[
        "biopython",
        "numpy",
        "scipy",
    ],
    python_requires=">=3.7",
)