from pathlib import Path

from setuptools import find_packages, setup


def read_readme() -> str:
	readme = Path(__file__).with_name("README.md")
	return readme.read_text(encoding="utf-8") if readme.exists() else ""


setup(
	name="inSTRbility",
	version="0.1.0",
	description="Installation package for inSTRbility",
	long_description=read_readme(),
	long_description_content_type="text/markdown",
	packages=find_packages(),
	include_package_data=True,
	install_requires=[],
)
