from setuptools import setup

setup(
    name="hvreglocator",
    version="0.1",
    py_modules=['hvreglocator'],
    entry_points={
        "console_scripts": [
            "hvreglocator=hvreglocator:main",
        ],
    },
    install_requires=[
        "biopython",
        "numpy",
        "scipy",
        "matplotlib"
    ],
    python_requires=">=3.7",
)