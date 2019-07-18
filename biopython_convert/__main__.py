#!/usr/bin/env python
import sys

from . import get_args, convert

def main():
    convert(*get_args(sys.argv[1:]))

if __name__ == "__main__":
    main()
