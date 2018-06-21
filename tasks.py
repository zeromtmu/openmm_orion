# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import os
import sys
import shutil
from invoke import task, run
from packaging.version import Version
import multiprocessing


from MDOrion import __version__  # noqa


@task
def flake8(ctx):
    run("flake8 MDOrion")


@task
def test(ctx):
    """
    run tests
    """
    clean_pyc(ctx)
    flake8(ctx)
    try:
        from xdist import __version__ as xver  # noqa
        os.system("python -m pytest -s -n {} --tb=native".format(multiprocessing.cpu_count()))
    except ImportError:
        os.system("python -m pytest -s --tb=native")


@task
def package(ctx):
    """
    Create package
    """
    clean(ctx)
    append_dev_version(ctx)
    run("python setup.py sdist --formats=gztar")
    run("python setup.py bdist_wheel")
    reset_dev_version(ctx)


@task
def upload(ctx):
    """
    Upload
    """
    current_ver = append_dev_version(ctx)
    package(ctx)
    twine_dist_args = "upload -r magpie ./dist/OpenEye-MDOrion-{}.tar.gz".format(current_ver)
    twine_wheel_args = "upload -r magpie ./dist/OpenEye_MDOrion-{}-*.whl".format(current_ver)
    if sys.platform == 'win32':
        run("twine.exe {}".format(twine_dist_args))
        run("twine.exe {}".format(twine_wheel_args))
    else:
        run("twine {}".format(twine_dist_args))
        run("twine {}".format(twine_wheel_args))


@task
def append_dev_version(ctx):
    """
    set dev version tag - git timestamp and sha
    """
    if "dev" not in __version__:
        return __version__

    print("Updating dev version with timestamp and git SHA")

    results = run("git log --pretty=format:'%cd+git%h' -n 1 --date=format:%Y%m%d%H%M%S --abbrev=8",
                  hide='both')
    dev_version_tag = results.stdout

    new_lines = []

    with open('MDOrion/__init__.py', 'r') as f:
        all_lines = f.readlines()
        for line in all_lines:
            if line.startswith("__version__"):
                vers = Version(line.split("=")[-1].strip().replace("'", "").replace('"', ''))
                new_lines.append('__version__ = "{}.dev{}"\n'.format(vers.base_version,
                                                                     dev_version_tag))
                dev_version = "{}.dev{}".format(vers.base_version, dev_version_tag)
            else:
                new_lines.append(line)

    with open('MDOrion/__init__.py', 'w') as o:
        o.writelines(new_lines)
        o.close()

    return dev_version


@task
def reset_dev_version(ctx):
    """
    reset/clean the dev version tag
    """
    if "dev" not in __version__:
        return

    print("Cleaning dev version with timestamp and git SHA")

    new_lines = []

    with open('MDOrion/__init__.py', 'r') as f:
        all_lines = f.readlines()
        for line in all_lines:
            if line.startswith("__version__"):
                vers = Version(line.split("=")[-1].strip().replace("'", "").replace('"', ''))
                new_lines.append('__version__ = "{}.dev"\n'.format(vers.base_version))
            else:
                new_lines.append(line)

    with open('MDOrion/__init__.py', 'w') as o:
        o.writelines(new_lines)
        o.close()


@task
def update_copyright(ctx):
    run("oecopyright --type example examples")
    run("oecopyright --type proprietary MDOrion")


@task
def docs(ctx):
    clean_docs(ctx)
    curdir = os.getcwd()
    os.chdir('docs')
    if sys.platform == 'win32':
        run("make.bat html")
    else:
        run("make html")
    os.chdir(curdir)


@task
def serve_docs(ctx):
    docs(ctx)
    curdir = os.getcwd()
    os.chdir('docs/build/html')
    run('python -m http.server')
    os.chdir(curdir)


@task
def clean(ctx):
    """
    Clean up doc and package builds
    """
    clean_pyc(ctx)
    clean_docs(ctx)
    clean_pycache(ctx)
    shutil.rmtree("dist", ignore_errors=True)
    shutil.rmtree("build", ignore_errors=True)
    egg_path = "{}.egg-info".format("MDOrion".replace("-", "_"))
    if os.path.isfile(egg_path):
        os.remove(egg_path)
    elif os.path.isdir(egg_path):
        shutil.rmtree(egg_path)


@task
def clean_pyc(ctx):
    """
    cleans out .pyc files
    """
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".pyc"):
                filename = os.path.join(root, file)
                if os.path.exists(filename):
                    os.unlink(filename)


@task
def clean_pycache(ctx):
    """
    cleans out __pycache__ dirs
    """
    for dirpath, dirs, files in os.walk(os.getcwd()):
        for dir in dirs:
            if dir == '__pycache__':
                del_dir = os.path.join(dirpath, dir)
                shutil.rmtree(del_dir)
@task
def clean_docs(ctx):
    doc_dir = "docs/build/html"
    _clean_out_dir(doc_dir)

    if os.path.isdir("docs/build/doctrees"):
        shutil.rmtree("docs/build/doctrees")


def _clean_out_dir(dir_path):
    if os.path.isdir(dir_path):
        for the_file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, the_file)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(e)
