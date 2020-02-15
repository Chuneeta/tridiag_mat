from setuptools import setup
import sys
import os
from tridiag_mat import version
import json

data = [version.git_origin, version.git_hash, version.git_description, version.git_branch]
with open(os.path.join('tridiag_mat', 'GIT_INFO'), 'w') as outfile:
    json.dump(data, outfile)

def package_files(package_dir, subdirectory):
    # walk the input package_dir/subdirectory
    # return a package_data list
    paths = []
    directory = os.path.join(package_dir, subdirectory)
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            path = path.replace(package_dir + '/', '')
            paths.append(os.path.join(path, filename))
    return paths

setup_args = {
    'name':         'tridiag_mat',
    'author':       'Chuneeta Nunhokee',
    'url':          'https://github.com/Chuneeta/beam_solver',
    'license':      'BSD',
    'version':      version.version,
    'description':  'HERA Primary Beam Estimator.',
    'packages':     ['tridiag_mat'],
    'package_dir':  {'tridiag_mat': 'tridiag_mat'},
    'install_requires': ['numpy>=1.14', 'matplotlib>=2.2'],
    'include_package_data': True,
    'zip_safe':     False,
}

if __name__ == '__main__':
    setup(*(), **setup_args)
# apply(setup, (), setup_args)
