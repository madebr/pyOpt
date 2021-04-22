#!/usr/bin/env python

import os


def configuration(parent_package, top_path):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('pyOpt', parent_package, top_path)

    # need: auto add_subpackage from source availability
    config.add_subpackage('pyALGENCAN')
    config.add_subpackage('pyALHSO')
    config.add_subpackage('pyALPSO')
    config.add_subpackage('pyCOBYLA')
    config.add_subpackage('pyCONMIN')
    config.add_subpackage('pyFILTERSD')
    config.add_subpackage('pyFSQP')
    config.add_subpackage('pyGCMMA')
    # config.add_subpackage('pyIPOPT')
    config.add_subpackage('pyKSOPT')
    config.add_subpackage('pyMIDACO')
    config.add_subpackage('pyMMA')
    config.add_subpackage('pyMMFD')
    config.add_subpackage('pyNLPQL')
    config.add_subpackage('pyNLPQLP')
    config.add_subpackage('pyNSGA2')
    config.add_subpackage('pyPSQP')
    config.add_subpackage('pySDPEN')
    config.add_subpackage('pySLSQP')
    config.add_subpackage('pySNOPT')
    config.add_subpackage('pySOLVOPT')

    for fn in 'LICENSE', 'README.md':
        config.add_data_files(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), fn))

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
