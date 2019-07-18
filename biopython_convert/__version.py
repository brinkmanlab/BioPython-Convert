from pbr.version import VersionInfo
_v = VersionInfo('mock').semantic_version()
__version__ = _v.release_string()

if __name__ == '__main__':
    print(__version__)
