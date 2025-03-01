# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run bridge unit tests for cclib."""

import sys
import unittest

sys.path.insert(1, "bridge")

if sys.version_info[1] >= 6:
    from .bridge.testpsi4 import *
    from .bridge.testpyscf import *
    from .bridge.testhorton import HortonTest
if sys.version_info[1] >= 5:
    from .bridge.testase import *
if sys.version_info[1] >= 4:
    from .bridge.testbiopython import *
from .bridge.testopenbabel import *
from .bridge.testpyquante import pyquante2Test

if __name__ == "__main__":
    unittest.main()
