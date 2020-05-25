#!/usr/bin/env python

import os
import sys
from numpy.distutils.command.build_ext import build_ext

if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

if sys.version_info[:2] < (2, 4):
    print(('pyOpt requires Python version 2.4 or later (%d.%d detected).' %sys.version_info[:2]))
    sys.exit(-1)


class build_opt(build_ext):
    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except:
            self.announce('*** WARNING: Building of optimizer "%s" '
            'failed: %s' %(ext.name, sys.exc_info()[1]))


def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration(None,parent_package,top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True,
    )
    
    config.add_subpackage('pyOpt')
    
    return config
    

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(
        name             = 'pyOpt',
        version          = '1.2.0',
        author           = 'Ruben E. Perez, Peter W. Jansen',
        author_email     = 'Ruben.Perez@rmc.ca; Peter.Jansen@rmc.ca',
        maintainer       = 'pyOpt Developers',
        maintainer_email = 'Ruben.Perez@rmc.ca; Peter.Jansen@rmc.ca',
        url              = 'http://pyopt.org/',
        download_url     = 'http://pyopt.org/',
        description      = 'Python package for formulating and solving nonlinear constrained optimization problems',
        long_description = 'pyOpt is a Python package for formulating and solving nonlinear constrained optimization problems',
        keywords         = 'optimization',
        license          = 'GNU LGPL',
        platforms        = ['Windows','Linux','Solaris','Mac OS-X','Unix'],
        classifiers      = ['Development Status :: 5 - Production/Stable',
                            'Environment :: Console',
                            'Intended Audience :: Science/Research',
                            'Intended Audience :: Developers',
                            'Intended Audience :: Education',
                            'License :: LGPL',
                            'Operating System :: Microsoft :: Windows',
                            'Operating System :: POSIX :: Linux',
                            'Operating System :: Unix',
                            'Operating System :: MacOS',
                            'Programming Language :: Python',
                            'Topic :: Scientific/Engineering',
                            'Topic :: Software Development',
                            'Topic :: Education'],
        configuration    = configuration,
        cmdclass = {"build_ext": build_opt},
    )
