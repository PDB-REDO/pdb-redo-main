"""
Helper functions for PDB_REDO python programs.
"""


from __future__ import division


import argparse
import bz2
import gzip
import logging
import math
import os
import sys


class DefaultHelpParser(argparse.ArgumentParser):
    def error(self, message):
        """Override error to show help by default."""
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


class Vector(object):
    """A simple class to represent vectors of three dimensions x, y, z."""

    def __init__(self, x, y, z):
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    def __repr__(self):
        return '<Vector x:{0:.3f}, y:{1:.3f}, z:{2:.3f}>'.format(self.x,
                                                                 self.y,
                                                                 self.z)

    def __str__(self):
        return '({0:.3f}, {1:.3f}, {2:.3f})'.format(self.x, self.y, self.z)

    def __add__(self, u):
        """Elementwise addition."""
        if isinstance(u, Vector):
            w = Vector(self.x + u.x, self.y + u.y, self.z + u.z)
        else:
            # scalar
            w = Vector(self.x + u, self.y + u, self.z + u)
        return w

    def __div__(self, u):
        """Elementwise division."""
        if isinstance(u, Vector):
            w = Vector(self.x / u.x, self.y / u.y, self.z / u.z)
        else:
            # scalar
            w = Vector(self.x - u, self.y - u, self.z - u)
        return w

    def __mul__(self, u):
        """Dot product."""
        return self.x*u.x + self.y*u.y + self.z*u.z

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)

    def __pow__(self, u):
        """Cross product.

        u x v = u2v3-u3v2, u3v1-u1v3, u1v2-v2u1
        """
        return Vector(self.y*u.z - self.z*u.y,
                      self.z*u.x - self.x*u.z,
                      self.x*u.y - self.y*u.x)

    def __sub__(self, u):
        """Elementwise subtraction."""
        if isinstance(u, Vector):
            w = Vector(self.x - u.x, self.y - u.y, self.z - u.z)
        else:
            # scalar
            w = Vector(self.x - u, self.y - u, self.z - u)
        return w

    def angle(self, u):
        """Angle in radians."""
        nv = self.l2_norm()
        nu = u.l2_norm()
        quot = (self * u) / (nv * nu)
        # Cos is defined between -1 and 1
        quot = min(quot, 1)
        quot = max(quot, -1)
        return math.acos(quot)

    def l1_norm(self):
        """L1 or Manhattan norm."""
        return abs(self.x) + abs(self.y) + abs(self.z)

    def l2_norm(self):
        """L2 or Euclidian norm."""
        return math.sqrt(self*self)

    def linf_norm(self):
        """L-infinity norm."""
        return max(abs(self.x), abs(self.y), abs(self.z))

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @classmethod
    def from_atom_record(cls, atom_record):
        x, y, z = atom_record[30:38], atom_record[38:46], atom_record[46:54]
        vector = cls(x, y, z)
        return vector


class VersionActionStdOut(argparse._VersionAction):
    def __call__(self, parser, namespace, values, option_string=None):
        """Override __call__ to print version to stdout."""
        version = self.version
        if version is None:
            version = parser.version
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), sys.stdout)
        parser.exit()


def is_valid_file(parser, arg, empty_allowed=False):
    """Return file name arg if file exists."""
    if not os.path.isfile(arg):
        parser.error('The file {} does not exist!'.format(arg))
    elif not empty_allowed and not os.stat(arg).st_size > 0:
        parser.error('The file {} is empty!'.format(arg))
    else:
        return arg


def read(pdb_file_path):
    """Return lines from uncompressed, .gz or .bz2 file."""
    if pdb_file_path.endswith('.gz'):
        with gzip.open(pdb_file_path) as fh:
            return fh.readlines()

    if pdb_file_path.endswith('.bz2'):
        with bz2.BZ2File(pdb_file_path) as fh:
            return fh.readlines()

    with open(pdb_file_path) as fh:
        return fh.readlines()


def setup_logger(log, verbose=False):
    log.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    if not verbose:
        handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    return handler


def write_line(line, path):
    """Write string line to uncompressed, .gz or .bz2 file."""
    if path.endswith('.gz'):
        with gzip.open(path, 'w') as fh:
            return fh.write(line)
    elif path.endswith('.bz2'):
        with bz2.BZ2File(path, 'w') as fh:
            return fh.write(line)
    else:
        with open(path, 'w') as fh:
            return fh.write(line)
