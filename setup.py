import setuptools
import os.path
from os import listdir
import re
from setuptools import setup
from pathlib import Path


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name='rfpy',
    version=find_version('rfpy', '__init__.py'),
    description='Python Module for Teleseismic Receiver Functions',
    author='Pascal Audet',
    author_email='pascal.audet@uottawa.ca',
    maintainer='Pascal Audet',
    maintainer_email='pascal.audet@uottawa.ca',
    url='https://github.com/paudetseis/RfPy',
    classifiers=[
         'Development Status :: 3 - Alpha',
         'License :: OSI Approved :: MIT License',
         'Programming Language :: Python :: 3.6',
         'Programming Language :: Python :: 3.7',
         'Programming Language :: Python :: 3.8',
         'Programming Language :: Python :: 3.9',
         'Programming Language :: Python :: 3.10'],
    install_requires=['numpy', 'obspy', 'stdb>=0.2.0'],
    python_requires='>=3.7',
    packages=setuptools.find_packages(),
    entry_points={'console_scripts': [
        'rfpy_calc=rfpy.scripts.rfpy_calc:main',
        'rfpy_recalc=rfpy.scripts.rfpy_recalc:main',
        'rfpy_plot=rfpy.scripts.rfpy_plot:main',
        'rfpy_harmonics=rfpy.scripts.rfpy_harmonics:main',
        'rfpy_hk=rfpy.scripts.rfpy_hk:main',
        'rfpy_ccp=rfpy.scripts.rfpy_ccp:main']})
