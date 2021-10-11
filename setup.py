import setuptools
import dubins


with open("README.md") as rm:
    long_description = rm.read()


setuptools.setup(
    name='dubins',
    version='0.0',
    author='Maksim Buzikov',
    author_email="shmax3@gmail.com",
    description="Toolkit for problems containing Dubins model",
    long_decription=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shmax3/Dubins",
    test_suite='dubins.tests',
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'scipy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)