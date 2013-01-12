#!/usr/bin/env python

import os,sys


def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('pyFILTERSD',parent_package,top_path)
    
    config.add_library('filtersd',
        sources=[os.path.join('source', '*.f')])
    config.add_extension('filtersd',
        sources=['source/f2py/filtersd.pyf'],
        libraries=['filtersd'])
    config.add_data_files('LICENSE','README')
    
    return config
    

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
    
