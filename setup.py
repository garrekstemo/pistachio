import setuptools

with open("README.md", 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name="Pistachio",
    version="0.8.0",
    author="Garrek Stemo",
    author_email="stemo.garrek_danneker.se3@ms.naist.jp",
    description="A transfer matrix algorithm and related packages",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/garrekds/pistachio",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Alpha",
        "License :: MIT Licence",
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.7'
)
