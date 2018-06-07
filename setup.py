#!/usr/bin/env python
import re, ast, os
from os.path import relpath, join
from setuptools import setup, find_packages

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

def get_reqs(reqs):
    return [str(ir.req) for ir in reqs]

try:
    install_reqs = get_reqs(parse_requirements("requirements.txt"))
except TypeError:
    try: #for pip >= 10
        from pip._internal.download import PipSession
    except ImportError: # for pip <= 9.0.3
        from pip.download import PipSesion
    install_reqs = get_reqs(
        parse_requirements("requirements.txt", session=PipSession())
    )


def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files

def get_version():
    _version_re = re.compile(r'__version__\s+=\s+(.*)')
    with open('OpenMMCubes/__init__.py', 'rb') as f:
        version = str(ast.literal_eval(_version_re.search(f.read().decode('utf-8')).group(1)))
        return version

setup(
    name="OpenMMCubes",
    version='0.2.2',
    packages=find_packages(include=['examples'], exclude=['tests*']),
    include_package_data=True,
    package_data={ 'examples' : find_package_data('examples/data', 'examples')},
    author="Christopher Bayly, Gaetano Calabro, Nathan M. Lim, John Chodera, ",
    author_email="bayly@eyesopen.com",
    description='Prepare complex for MD with OpenMM',
    install_requires=install_reqs,
    license='Other/Proprietary License',
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
    ]
)
