import subprocess
import datetime
import os
import importlib
import hashlib


def _develop_packages(conda_list):
    dev_packages = []
    not_comment = lambda s: s[0][0] != '#'
    has_channel = lambda s: len(s) == 4
    is_develop = lambda s: s[3] == '<develop>'
    is_dev_package = (lambda s: has_channel(s) and not_comment(s)
                      and is_develop(s))
    for line in conda_list:
        splitted = line.split()
        if is_dev_package(splitted):
            dev_packages.append(splitted[0])
    return dev_packages

def output_file_hashes(output_files):
    pass

def get_conda_list():
    results = subprocess.run(["conda", "list"], capture_output=True)
    conda_list = results.stdout.decode('utf-8').split('\n')
    return conda_list

def get_development_package_versions(conda_list):
    dev_packages = {}
    develop_packages = _develop_packages(conda_list)
    for package in develop_packages:
        package_cleaned = package.translate(str.maketrans('-', '_'))
        pkg = importlib.import_module(package_cleaned)
        try:
            pkg_version = pkg.version.full_version
        except AttributeError:
            pkg_version = 'Unknown'

        dev_packages[package] = pkg_version

    return dev_packages

def dev_package_versions_string(dev_packages):
    version_strings = [
        package + ": " + version
        for package, version in dev_packages.items()
    ]
    strings = ["Full version for <develop> packages:"] + version_strings
    return strings

def reproducibility_string():
    timestamp = [str(datetime.datetime.now())]
    conda_list = get_conda_list()
    dev_packages = get_development_package_versions(conda_list)
    dev_package_list = dev_package_versions_string(dev_packages)
    return ('\n'.join(timestamp + conda_list + dev_package_list))

if __name__ == "__main__":
    print(reproducibility_string())
