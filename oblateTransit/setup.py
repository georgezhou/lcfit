#! /usr/bin/env python

from distutils.core import *
from distutils import sysconfig
from distutils.extension import Extension

import numpy

numpy_include = numpy.get_include()
_oblateness = Extension("_oblateness",
        ["oblateness_wrap.cxx",
            "oblateness.cc"],
        include_dirs = [numpy_include],
        )

setup(name="Oblate",
        description = "model for an oblate planet",
        author = "Wei Zhu, Xu Huang",
        author_email = "",
        url = "",
        version = "0.0.0",
        py_modules = ["oblateness"],
        ext_modules = [_oblateness])
