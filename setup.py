from setuptools import setup

setup(
    name="hvrlocator",
    version="0.2",
    py_modules=['hvrlocator'],
    entry_points={
        "console_scripts": [
            "hvrlocator=hvrlocator:main",
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
