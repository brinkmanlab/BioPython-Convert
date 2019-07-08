#!/usr/bin/env python
import sys

from . import get_args, convert

convert(*get_args(sys.argv[1:]))