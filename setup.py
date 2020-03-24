import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="opyrators",
    version="0.0.1",
    author="Evert van Nieuwenburg",
    author_email="evert.v.nieuwenburg@gmail.com",
    description="A minimal package to compute commutators of quantum operators"
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/everthemore/opyrators",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
)
