#!/usr/bin/env python

from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('stats', parent_package, top_path)

    config.add_data_dir('tests')

    # add lmoments module
    config.add_extension('_lmoments', sources=['lmoments.f', 'lmomextras.f', 'regionalization.f'], )
    
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
