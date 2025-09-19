from setuptools import setup, find_packages

setup(
    name="decid_package",
    version="0.1",
    packages=find_packages("scripts"),
    package_dir={"": "scripts"},
)
