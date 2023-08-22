#!/usr/bin/env python

import os
import sys

from numpy.distutils.command.build_ext import build_ext

if sys.version_info[:2] < (3, 7):
    raise RuntimeError('pyOpt requires Python version 3.7 or later ({:d}.{:d} detected).'.format(
        sys.version_info[:2]))
    sys.exit(-1)


class build_opt(build_ext):
    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except:
            self.announce(f'*** WARNING: Building of optimizer {ext.name} failed: {sys.exc_info()[1]}')


def configuration(parent_package='', top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True,
    )

    config.add_subpackage('pyOpt')

    return config

extras={
    "test": ["pytest"],
}
extras['dev'] = extras['test']

extras["all"] = sum(extras.values(), [])


if __name__ == '__main__':
    from numpy.distutils.core import setup

    # from distutils.command.sdist import sdist
    setup(
        configuration=configuration,
        cmdclass = {"build_ext": build_opt,
                    # 'sdist': sdist  # TODO: why does this not work? Why need manifest to include pyOpt?
                    },
        extras_require=extras

    )
