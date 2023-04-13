from setuptools import Command, find_packages, setup

__lib_name__ = "PROST"
__lib_version__ = "1.1.2"
__description__ = "PROST: A quantitative pattern recognition framework for spatial transcriptomics"
__url__ = "https://github.com/Tang-Lab-super/PROST"
__author__ = "Yuchen Liang"
__author_email__ = "liangych35@mail2.sysu.edu.cn"
__license__ = "MIT"
__keywords__ = ["spatial transcriptomics", "PROST index", "PROST Neural Network"]
__requires__ = ["requests",]

with open("README.rst", "r", encoding="utf-8") as f:
    __long_description__ = f.read()

setup(
    name = __lib_name__,
    version = __lib_version__,
    description = __description__,
    url = __url__,
    author = __author__,
    author_email = __author_email__,
    license = __license__,
    packages = ['PROST'],
    install_requires = __requires__,
    zip_safe = False,
    include_package_data = True,
    long_description = __long_description__
)
