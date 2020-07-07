import setuptools

with open("README.md", 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name="Pistachio",
    version="0.7.1",
    author="Garrek Stemo",
    author_email="stemo.garrek_danneker.se3@ms.naist.jp",
    description="A transfer matrix algorithm and related packages",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/garrekds/pistachio",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Alpha",
        "License :: MIT Licence",
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.7'
)
