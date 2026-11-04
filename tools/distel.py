#!/usr/bin/python

"""
  Version 1.03 2026-11-04

  Calculate  distance and angle restraint violation statistics.
  List rmsZ and outliers.

  Written by Robbie Joosten, Bart van Beusekom, Wouter Touw, Daniel Alvarez Salmoral
  E-mail:    r.joosten@nki.nl

  If you publish results (directly or indirectly) obtained by using
  distel: Please, refer to (one of) these references:
  - Robbie P. Joosten, Gert Vriend: "PDB improvement starts with data
    deposition" Science, 317, p. 195-196 (2007)
  - Robbie P. Joosten, Thomas Womack, Gert Vriend and Gerard Bricogne:
    "Re-refinement fromdeposited X-ray data can deliver improved models
    for most PDB entries"  Acta Cryst. D65, p. 176-185 (2009)
  - Robbie P. Joosten, Jean Salzemann, Vincent Bloch, Heinz Stockinger,
    Ann-Charlott Berglund, Christophe Blanchet, Erik Bongcam-Rudloff,
    Christophe Combet, Ana L. Da Costa, Gilbert Deleage, Matteo
    Diarena, Roberto Fabbretti, Geraldine Fettahi, Volker Flegel,
    Andreas Gisel, Vinod Kasam, Timo Kervinen, Eija Korpelainen, Kimmo
    Mattila, Marco Pagni, Matthieu Reichstadt, Vincent Breton, Ian J.
    Tickle, Gert Vriend: "PDB_REDO: automated re-refinement of X-ray
    structure models in the PDB" J. Appl. Cryst., 42, p. 376-384 (2009)
  - Robbie P. Joosten, Tim A.H. te Beek, Elmar Krieger, Maarten
    Hekkelman, Rob W.W. Hooft, Reinhard Schneider, Chris Sander, Gert
    Vriend: "A series of PDB related databases for everyday needs"
    Nucl. Acids Res., 39, p. D411-D419 (2011)
  - Robbie P. Joosten, Krista Joosten, Serge X. Cohen, Gert Vriend,
    Anastassis Perrakis: "Automatic rebuilding and optimization of
    crystallographic structures in the Protein Data Bank"
    Bioinformatics, 27, p. 3392-3398 (2011)
  - Robbie P. Joosten, Krista Joosten, Garib N. Murshudov, Anastassis
    Perrakis: "PDB_REDO: constructive validation, more than just
    looking for errors" Acta Cryst. D68, p. 484-496 (2012)

  Change log
  Version 1.03
  - For Na and Mg sites the number of ligands, the mean coordination length and the BV score is given
  Version 1.02
  - Apply the symmetry operator from $CLIBD/symop.lib to angle restraints
    with 'symm Y' instead of ignoring them.
  Version 1.01
  - Add angle restraint statistics (rmsZ and outliers).
  - Report and ignore angle restraints that request symmetry ('symm Y').
  Version 1.00
  - Major rewrite using Big Pickle. 
  - Massive reduction in running time from code optimisation.
  - Additional speedup by (small) algorithmic change in the jackknifing.
  - Consolidated code from helpers.py.
  - Commented out orphaned code.
  Version 0.09
  - Bugfix in how the order of the atom_site columns is established.
  Version 0.08
  - Stripped program to bare minimum.
  Version 0.07
  - Now works with mmCIF files only.
  Version 0.06: 
  - Calculates Jackknife stdev estimate
  Version 0.05:
  - Add torsion angle restraint violations
  - Fix distance restraint rmsZ
  Version 0.04:
  - Fix restraint matching for rebuilt models
  Version 0.03:
  - Add logging and update style
  - Allow compressed PDB files
  Version 0.02:
  - If there are multiple targets for the same interatomic distance,
    only the smalles violation is taken into account.
  Version 0.01:
  - Writing the YASARA macro is now optional.
  Version 0.00:
  - A first attempt
"""


from __future__ import division

import argparse
import bz2
import gzip
import itertools
import logging
import math
import os
import re
import sys

import numpy as np

VERSION = 1.03

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

    def __mul__(self, u):
        """Dot product."""
        return self.x*u.x + self.y*u.y + self.z*u.z

    def __sub__(self, u):
        """Elementwise subtraction."""
        return Vector(self.x - u.x, self.y - u.y, self.z - u.z)

    def l2_norm(self):
        """L2 or Euclidian norm."""
        return math.sqrt(self*self)

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
    def from_atom_record(cls, atom_record, atom_site_columns):
        column_dict = {col: i for i, col in enumerate(atom_site_columns)}
        x, y, z = atom_record.split()[column_dict['Cartn_x']], \
                  atom_record.split()[column_dict['Cartn_y']], \
                  atom_record.split()[column_dict['Cartn_z']]
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


def read(mmcif_file_path):
    """Return lines from uncompressed, .gz or .bz2 file."""
    if mmcif_file_path.endswith('.gz'):
        with gzip.open(mmcif_file_path, 'rt',  encoding='utf-8') as fh:
            return fh.readlines()

    if mmcif_file_path.endswith('.bz2'):
        with bz2.open(mmcif_file_path, 'rt',  encoding='utf-8') as fh:
            return fh.readlines()

    with open(mmcif_file_path, 'rt',  encoding='utf-8') as fh:
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


LOG = logging.getLogger(__name__)


def _main():
    """Compare restraints and structure.

    Restraint files must have the external restraints format (REFMAC).
    Distance and angle restraints have been implemented.
    """
    fmt_version = 'distel (version {})'
    descr = '{} - {}'.format(fmt_version.format(VERSION),
                             'Analyze restraints.')
    parser = DefaultHelpParser(description=descr)
    parser.add_argument('-v', '--verbose', help='Verbose mode',
                        action='store_true')
    parser.add_argument('--symop', help='Path to symop.lib '
                        '(default: $CLIBD/symop.lib)', default=None)
    parser.add_argument('--version', help='Print version',
                        action=VersionActionStdOut,
                        version=fmt_version.format(VERSION))
    parser.set_defaults(verbose=False, print_version=False)
    parser.add_argument('mmcif_file_path',
                        help='Path to input mmCIF file (.gz and .bz2 supported)',
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('rest_file_path',
                        help='Path to input restraint file '
                             '(.gz and .bz2 supported)',
                        type=lambda x: is_valid_file(parser, x,
                                                     empty_allowed=True))
    args = parser.parse_args()

    setup_logger(LOG, args.verbose)

    space_group = read_space_group(args.mmcif_file_path)

    restraint_pairs = read_restraints(args.rest_file_path)
    distance_restraints, angle_restraints = restraint_pairs
    mmcif_atoms, atom_site_columns = read_mmcif(args.mmcif_file_path)
    distance_matches = match_restraints_to_atoms(distance_restraints,
                                                 mmcif_atoms,
                                                 atom_site_columns)

    calc_distance_deviations(distance_matches, atom_site_columns)

    if angle_restraints:
        angle_matches = match_angle_restraints_to_atoms(angle_restraints,
                                                        mmcif_atoms,
                                                        atom_site_columns)
        symmetry = None
        if any(any(restraint[4]) for restraint in angle_restraints):
            symop_file = args.symop or os.path.join(
                os.environ.get('CLIBD', ''), 'symop.lib')
            try:
                symmetry = read_symmetry(args.mmcif_file_path, symop_file)
            except ValueError as error:
                sys.exit('error: {}'.format(error))
        calc_angle_deviations(angle_matches, atom_site_columns, symmetry)
        calc_metal_ligand_summary(angle_matches, atom_site_columns, symmetry)
    
    
# Unused: superseded by the closed-form jackknife_stdev below.
# def jackknife(x, func):
#     """Jackknife estimate of the estimator func"""
#     n = len(x)
#     idx = np.arange(n)
#     return sum(func(x[idx!=i]) for i in range(n))/float(n)

def jackknife_stdev(x):
    """Jackknife estimate of the variance of the rms estimator.

    Closed-form leave-one-out computation, O(n) instead of O(n^2).
    """
    n = len(x)
    if n < 2:
        return np.nan
    x = np.asarray(x, dtype=float)
    sum_sq = float((x * x).sum())
    # rms of x with element i left out, for every i
    leave_one_out = np.sqrt((sum_sq - x * x) / (n - 1))
    j_est = leave_one_out.mean()
    return np.sqrt((n - 1) / float(n) * ((leave_one_out - j_est) ** 2).sum())

# Unused: the leave-one-out rms is computed directly in jackknife_stdev.
# def rms(x):
#     return np.sqrt(np.mean(x**2))
    

def calc_distance(v1, v2):
    """Return Euclidian distance between two vectors."""
    dif = v1 - v2
    return dif.l2_norm()


def format_atom(atom_record, atom_site_columns):
    """Format atom identity as auth_asym_id auth_seq_id ins_code auth_atom_id label_alt_id."""
    column_dict = {col: i for i, col in enumerate(atom_site_columns)}
    words = atom_record.split()
    auth_asym_id = words[column_dict['auth_asym_id']]
    auth_seq_id = words[column_dict['auth_seq_id']]
    inscode = words[column_dict['pdbx_PDB_ins_code']]
    if inscode == '?':
        inscode = '.'
    atom_name = words[column_dict['auth_atom_id']]
    alt_code = words[column_dict['label_alt_id']]
    if alt_code == '.':
        alt_code = '.'
    return '{0:s} {1:s} {2:s} {3:s} {4:s}'.format(auth_asym_id, auth_seq_id,
                                                  inscode, atom_name,
                                                  alt_code)


def calc_distance_deviations(matches, atom_site_columns):
    """Calculate rmsZ and outlier statistics for distance restraints.

    Return None if there are no matches or all restraints are ignored.
    """
    if len(matches) == 0:
        LOG.info('Distance restraint rmsZ: NA NA')
        return None

    # Calculate the deviation for every match once. If there are multiple
    # restraints for the same atom pair, only the smallest violation is
    # taken into account (ties are all kept).
    deviations = []
    pair_indices = {}
    for i, match in enumerate(matches):
        target = match[2][0]
        sigma = match[2][1]
        dist = calc_distance(Vector.from_atom_record(match[0], atom_site_columns),
                             Vector.from_atom_record(match[1], atom_site_columns))
        deviations.append((dist, (dist - target) / sigma))
        pair_indices.setdefault((match[0], match[1]), []).append(i)

    keep = [True] * len(matches)
    for indices in pair_indices.values():
        if len(indices) > 1:
            best = min(abs(deviations[i][1]) for i in indices)
            for i in indices:
                if abs(deviations[i][1]) != best:
                    keep[i] = False

    total_z_sq, num_restraints_ignored = 0.0, 0
    # Unused: remnants of the optional YASARA macro output.
    # yasara_lines = []
    lzlength = []
    for i, match in enumerate(matches):
        dist, z = deviations[i]
        if keep[i]:
            lzlength.append(z)
            total_z_sq += z**2
            target = match[2][0]
            sigma = match[2][1]
            msg = '{0:s} -- {1:s} | {2:6.2f} {3:6.2f} {4:6.2f} ' \
                  '{5:6.2f}'.format(format_atom(match[0], atom_site_columns),
                                    format_atom(match[1], atom_site_columns),
                                    target, sigma, dist, z)
            if abs(z) > 4:
                LOG.info('OUTLIER: %s', msg)
            else:
                LOG.debug('         %s', msg)
        else:
            num_restraints_ignored += 1

    azlength = np.array(lzlength)
    #jrmsz = jackknife(azlength, rms)
    stdev = jackknife_stdev(azlength)
    #print (jrmsz, stdev)
    if len(matches) == num_restraints_ignored:
        LOG.info('All restraints ignored')
        LOG.info('Distance restraint rmsZ: NA NA')
        return None
    rmsz = math.sqrt(total_z_sq / (len(matches) - num_restraints_ignored))
    LOG.info('Distance restraint rmsZ: %6.3f %6.3f', rmsz, stdev)
    # return yasara_lines


def calc_angle(v1, v2, v3):
    """Return the angle in degrees between v1-v2 and v3-v2."""
    u = v1 - v2
    w = v3 - v2
    quot = (u * w) / (u.l2_norm() * w.l2_norm())
    # Cos is defined between -1 and 1
    quot = max(-1.0, min(1.0, quot))
    return math.degrees(math.acos(quot))


def calc_angle_deviations(matches, atom_site_columns, symmetry=None):
    """Calculate rmsZ and outlier statistics for angle restraints.

    Angles are reported in degrees. Atoms with 'symm Y' are placed with the
    symmetry data (see read_symmetry). Return None if there are no matches
    or all restraints are ignored.
    """
    if len(matches) == 0:
        LOG.info('Angle restraint rmsZ: NA NA')
        return None

    # Calculate the deviation for every match once. If there are multiple
    # restraints for the same atom triple, only the smallest violation is
    # taken into account (ties are all kept).
    deviations = []
    triple_indices = {}
    for i, match in enumerate(matches):
        target = match[3][0]
        sigma = match[3][1]
        ang = calc_angle(*place_angle_atoms(match, atom_site_columns, symmetry))
        deviations.append((ang, (ang - target) / sigma))
        triple_indices.setdefault((match[0], match[1], match[2]), []).append(i)

    keep = [True] * len(matches)
    for indices in triple_indices.values():
        if len(indices) > 1:
            best = min(abs(deviations[i][1]) for i in indices)
            for i in indices:
                if abs(deviations[i][1]) != best:
                    keep[i] = False

    total_z_sq, num_restraints_ignored = 0.0, 0
    lzlength = []
    for i, match in enumerate(matches):
        ang, z = deviations[i]
        if keep[i]:
            lzlength.append(z)
            total_z_sq += z**2
            target = match[3][0]
            sigma = match[3][1]
            msg = '{0:s} -- {1:s} -- {2:s} | {3:6.2f} {4:6.2f} ' \
                  '{5:6.2f} {6:6.2f}'.format(format_atom(match[0], atom_site_columns),
                                             format_atom(match[1], atom_site_columns),
                                             format_atom(match[2], atom_site_columns),
                                             target, sigma, ang, z)
            if abs(z) > 4:
                LOG.info('OUTLIER: %s', msg)
            else:
                LOG.debug('         %s', msg)
        else:
            num_restraints_ignored += 1

    azlength = np.array(lzlength)
    stdev = jackknife_stdev(azlength)
    if len(matches) == num_restraints_ignored:
        LOG.info('All angle restraints ignored')
        LOG.info('Angle restraint rmsZ: NA NA')
        return None
    rmsz = math.sqrt(total_z_sq / (len(matches) - num_restraints_ignored))
    LOG.info('Angle restraint rmsZ: %6.3f %6.3f', rmsz, stdev)


def place_angle_atoms(match, atom_site_columns, symmetry):
    """Return the three Vectors of an angle match, applying symmetry.

    Each atom with 'symm Y' is moved to the symmetry image closest to a
    reference atom. The reference is the middle atom, or the first atom
    without 'symm Y' if the middle atom has it.
    """
    flags = match[4]
    vectors = [Vector.from_atom_record(match[i], atom_site_columns)
               for i in range(3)]
    if not any(flags):
        return vectors
    ref = 1 if not flags[1] else next((i for i in (0, 2) if not flags[i]), 1)
    operators, matrix = symmetry
    inverse = np.linalg.inv(matrix)
    ref_frac = inverse.dot(vector_xyz(vectors[ref]))
    for i in range(3):
        if flags[i] and i != ref:
            atom_frac = inverse.dot(vector_xyz(vectors[i]))
            image = nearest_symm_image(ref_frac, atom_frac, operators, matrix)[3]
            vectors[i] = Vector(*matrix.dot(image))
    return vectors


# Bond valence parameters (Ro, B) for bond valence = exp((Ro - R) / B).
# Source: bvparm2020.cif (I. D. Brown, accumulated table of bond valence
# parameters, 2020), _valence_param_ref_id 'bt': Yee, Suescun & Rabuffetti
# (2018) J. Solid State Chem. 270, 242-246. Rows: Mg 2 O -2 1.611 0.37 and
# Na 1 O -2 1.744 0.37.
BOND_VALENCE_PARAMS = {
    ('Mg', 'O'): (1.611, 0.37),
    ('Mg', 'N'): (1.87, 0.37),
    ('Na', 'O'): (1.744, 0.37),
    ('Na', 'N'): (1.88, 0.37),
}


def calc_metal_ligand_summary(matches, atom_site_columns, symmetry=None):
    """Log the number of ligands, the mean distance and the bond valence sum.

    The ligands of a metal are the other atoms of the angle restraints the
    metal occurs in. Ligands with 'symm Y' use the closest symmetry image.
    The bond valence sum uses BOND_VALENCE_PARAMS for each metal-ligand
    element pair and is NA if no parameter is available.
    """
    column_dict = {col: i for i, col in enumerate(atom_site_columns)}
    sites = {}
    for match in matches:
        flags = match[4]
        for i in range(3):
            element = match[i].split()[column_dict['type_symbol']]
            if element.upper() not in ('NA', 'MG') or flags[i]:
                continue
            _, ligands = sites.setdefault(format_atom(match[i],
                                                      atom_site_columns),
                                          (match[i], {}))
            for j in range(3):
                if j != i:
                    ligands[format_atom(match[j], atom_site_columns)] = \
                        (match[j], flags[j])
    if not sites:
        return

    LOG.info('NA and MG summary: atom ligands mean_distance bond_valence')
    missing = set()
    for metal_line, ligands in sites.values():
        metal_xyz = Vector.from_atom_record(metal_line, atom_site_columns)
        metal_element = metal_line.split()[column_dict['type_symbol']]
        distances, valence, has_param = [], 0.0, False
        for ligand_line, symm in ligands.values():
            ligand_xyz = Vector.from_atom_record(ligand_line, atom_site_columns)
            if symm:
                operators, matrix = symmetry
                inverse = np.linalg.inv(matrix)
                dist = nearest_symm_image(
                    inverse.dot(vector_xyz(metal_xyz)),
                    inverse.dot(vector_xyz(ligand_xyz)),
                    operators, matrix)[0]
            else:
                dist = calc_distance(metal_xyz, ligand_xyz)
            distances.append(dist)
            ligand_element = ligand_line.split()[column_dict['type_symbol']]
            pair = (metal_element.capitalize(), ligand_element.capitalize())
            if pair in BOND_VALENCE_PARAMS:
                ro, b = BOND_VALENCE_PARAMS[pair]
                valence += math.exp((ro - dist) / b)
                has_param = True
            else:
                missing.add(pair)
        bond_valence = '{:6.3f}'.format(valence) if has_param else 'NA'
        LOG.info('%s %d %6.3f %s', format_atom(metal_line, atom_site_columns),
                 len(distances), sum(distances) / len(distances), bond_valence)
    for metal, ligand in sorted(missing):
        sys.stderr.write('No bond valence parameter for {}-{} (ref bt)\n'.format(
            metal, ligand))


def index_atoms(mmcif_atoms, atom_site_columns):
    """Return a dict mapping atom identity to ATOM/HETATM line.

    Building this index once turns restraint matching from
    O(restraints * atoms) into O(atoms + restraints).
    """
    column_dict = {col: i for i, col in enumerate(atom_site_columns)}
    atom_index = {}
    for atom_line in mmcif_atoms:
        atom_index[parse_atom_key(atom_line.split(), column_dict)] = atom_line
    return atom_index


def match_restraints_to_atoms(restraint_tuples, mmcif_atoms, atom_site_columns):
    """Determine ATOM/HETATM lines that match restraint atom tuples.

    Return a list of (atom_line_1, atom_line_2, (target, sigma)) tuples.
    """
    atom_index = index_atoms(mmcif_atoms, atom_site_columns)

    two_matches = []
    for restraint in restraint_tuples:
        match1 = atom_index.get(restraint[0])
        match2 = atom_index.get(restraint[1])
        if match1 and match2 and match1 is not match2:
            two_matches.append((match1, match2, restraint[2]))

    LOG.debug('%d restraints matched with two atoms from model data',
              len(two_matches))
    return two_matches


def match_angle_restraints_to_atoms(restraint_tuples, mmcif_atoms,
                                    atom_site_columns):
    """Determine ATOM/HETATM lines that match angle restraint atom tuples.

    Return a list of (atom_line_1, atom_line_2, atom_line_3,
    (target, sigma), symm) tuples, where symm holds the per-atom symmetry
    flags.
    """
    atom_index = index_atoms(mmcif_atoms, atom_site_columns)

    three_matches = []
    for restraint in restraint_tuples:
        match1 = atom_index.get(restraint[0])
        match2 = atom_index.get(restraint[1])
        match3 = atom_index.get(restraint[2])
        if match1 and match2 and match3 and \
                match1 is not match2 and match1 is not match3 and \
                match2 is not match3:
            three_matches.append((match1, match2, match3, restraint[3],
                                  restraint[4]))

    LOG.debug('%d restraints matched with three atoms from model data',
              len(three_matches))
    return three_matches

def read_mmcif(mmcif_file):
    """Return a list of ATOM and HETATM lines from model file."""
    lines = read(mmcif_file)
    atom_site_columns = []
    atoms = []
    capture = False
    previous_line_was_loop = False
    
    for line in lines:
        line = line.strip()
        if previous_line_was_loop and line.startswith('_atom_site.'):
            capture = True
            atom_site_columns.append(line.split('.')[1])  # Add the current line
        elif capture:
            if line.startswith('_atom_site.'):
                atom_site_columns.append(line.split('.')[1])  # Continue adding
            else:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atoms.append(line)
        
        # Check if current line is 'loop_' for the next iteration
        previous_line_was_loop = line == 'loop_'
        
    
    LOG.debug('%d atoms read.', len(atoms))
    return atoms, atom_site_columns


def read_space_group(mmcif_file):
    """Return the space group name from _symmetry.space_group_name_H-M."""
    lines = read(mmcif_file)
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('_symmetry.space_group_name_H-M'):
            # Value may be on the same line, possibly quoted
            parts = stripped.split(None, 1)
            if len(parts) == 2:
                value = parts[1].strip()
            else:
                continue
            # Remove surrounding quotes
            if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
                value = value[1:-1]
            return value
    return None


CELL_KEYS = ['_cell.length_a', '_cell.length_b', '_cell.length_c',
             '_cell.angle_alpha', '_cell.angle_beta', '_cell.angle_gamma']


def vector_xyz(vector):
    """Return the Cartesian coordinates of a Vector as a numpy array."""
    return np.array([vector.x, vector.y, vector.z])


def read_cell(mmcif_file):
    """Return the unit cell parameters (a, b, c, alpha, beta, gamma)."""
    values = {}
    for line in read(mmcif_file):
        parts = line.split()
        if len(parts) == 2 and parts[0] in CELL_KEYS:
            values[parts[0]] = float(parts[1].split('(')[0])
    return [values[key] for key in CELL_KEYS]


def orthogonalization_matrix(cell):
    """Return the matrix converting fractional to Cartesian coordinates."""
    a, b, c, alpha, beta, gamma = cell
    alpha, beta, gamma = (math.radians(angle) for angle in (alpha, beta, gamma))
    cos_a, cos_b, cos_g, sin_g = \
        math.cos(alpha), math.cos(beta), math.cos(gamma), math.sin(gamma)
    cy = (cos_a - cos_b * cos_g) / sin_g
    cz = math.sqrt(1.0 - cos_b ** 2 - cy ** 2)
    return np.array([[a, b * cos_g, c * cos_b],
                     [0.0, b * sin_g, c * cy],
                     [0.0, 0.0, c * cz]])


def normalize_space_group(name):
    """Return a space group name without spaces and in upper case."""
    return name.replace(' ', '').upper()


def read_symop_lib(symop_file):
    """Return a dict mapping normalized space group names to operators.

    Each operator is a string like 'X,Y,Z' as listed in symop.lib.
    """
    groups = {}
    current = None
    with open(symop_file, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if ',' in line:
                if current is not None:
                    current.append(line.strip())
                continue
            names = re.findall(r"'([^']*)'", line)
            current = [] if names else None
            for name in names:
                groups.setdefault(normalize_space_group(name), current)
    return groups


def parse_symop(operator):
    """Return (rotation, translation) for an operator string like '-Y,X-Y,Z'."""
    rotation, translation = [], []
    for component in operator.replace(' ', '').split(','):
        row, offset = [0.0, 0.0, 0.0], 0.0
        for sign, term in re.findall(r'([+-]?)(\d+/\d+|\d+|[XYZ])', component):
            factor = -1.0 if sign == '-' else 1.0
            if term in ('X', 'Y', 'Z'):
                row['XYZ'.index(term)] += factor
            elif '/' in term:
                numerator, denominator = term.split('/')
                offset += factor * float(numerator) / float(denominator)
            else:
                offset += factor * float(term)
        rotation.append(row)
        translation.append(offset)
    return np.array(rotation), np.array(translation)


def nearest_symm_image(reference_frac, atom_frac, operators, matrix):
    """Return the symmetry image of an atom closest to a reference atom.

    Coordinates are fractional. Return (distance, operator, shift, image),
    where operator is the 1-based index of the operator in the list, shift
    the lattice translation applied and image the fractional coordinates.
    """
    best = None
    for index, (rotation, translation) in enumerate(operators, 1):
        symm_frac = rotation.dot(atom_frac) + translation
        base = -np.round(symm_frac - reference_frac)
        for offset in itertools.product((-1, 0, 1), repeat=3):
            shift = base + np.array(offset)
            image = symm_frac + shift
            dist = float(np.linalg.norm(matrix.dot(image - reference_frac)))
            if best is None or dist < best[0]:
                best = (dist, index, tuple(int(s) for s in shift), image)
    return best


def read_symmetry(mmcif_file, symop_file):
    """Return (operators, matrix) for the space group of an mmCIF file.

    operators are (rotation, translation) pairs from symop_file and matrix
    converts fractional to Cartesian coordinates.
    """
    space_group = read_space_group(mmcif_file)
    groups = read_symop_lib(symop_file)
    key = normalize_space_group(space_group) if space_group else None
    if key not in groups:
        raise ValueError('no symmetry operators found for space group {} '
                         'in {}'.format(space_group, symop_file))
    operators = [parse_symop(op) for op in groups[key]]
    return operators, orthogonalization_matrix(read_cell(mmcif_file))


# Unused: distance restraints are parsed by parse_distance_restraint.
# def parse_complete_atom(atom_words):
#     """Parse the list atom_words in external restraint format.
#
#     Return a tuple (residue_id, alt_loc, atom_name)
#     Return a ValueError if the list is not of the expected length or if the
#     a string cannot be cast to the expected type.
#     """
#     exp = 'chain <chain> resi <resi> inse <inse> atom <atom> alte <alte>'
#     if not len(atom_words) == 10:
#         raise ValueError('Expected atom format: {}'.format(exp))
#     ins, alte = atom_words[5], atom_words[9]
#     ins = ' ' if ins == '.' else ins
#     alte = ' ' if alte == '.' else alte
#     return (atom_words[1],int(atom_words[3]),ins), alte, atom_words[7]


def parse_distance_restraint(restraint_line):
    """Parse the atoms, target and sigma from distance restraint string.

    Return a tuple of (id_atom_1, id_atom_2, (target, sigma)) where
    id_atom is a string tuple of (residue_id, alt_loc, atom_name).
    Return None if the distance restraint cannot be parsed.
    """
    if not restraint_line.startswith('exte dist'):
        return None
    is_next_chain, is_next_resnum, is_next_atom, is_next_ins = [False]*4
    is_next_altcode, is_next_target, is_next_sigma = [False]*3
    ch1, resn1, altc1, ins1, atm1, ch2, resn2, altc2, ins2, atm2 = [None]*10
    ch, resn, atm, targ, sigma = [None]*5
    altc, ins = ' ', ' '
    for word in restraint_line.split():
        if word == 'chain':
            is_next_chain = True
        elif word == 'resi' or word == 'residue':
            is_next_resnum = True
        elif word == 'atom':
            is_next_atom = True
        elif word == 'ins':
            is_next_ins = True
        elif word == 'alte':
            is_next_altcode = True
        elif word == 'value':
            is_next_target = True
        elif word == 'sigma':
            is_next_sigma = True
        elif word == 'second':
            if not ch or not resn or not atm:
                return None
            # set data of first atom if all necessary data is present in the
            # restraint
            ch1 = ch
            resn1 = resn
            altc1 = altc
            ins1 = ins
            atm1 = atm
            # reset params before reading 2nd atom
            ch = None
            resn = None
            altc = ' '
            ins = ' '
            atm = None
        elif is_next_chain:
            is_next_chain = False
            ch = word
        elif is_next_resnum:
            is_next_resnum = False
            resn = int(word)
        elif is_next_atom:
            is_next_atom = False
            atm = word
        elif is_next_ins:
            is_next_ins = False
            if not word == '.':
                ins = word
        elif is_next_altcode:
            is_next_altcode = False
            if not word == '.':
                altc = word
        elif is_next_target:
            is_next_target = False
            targ = float(word)
        elif is_next_sigma:
            is_next_sigma = False
            sigma = float(word)
    # check if second atom has been set properly
    if not ch or not resn or not atm or not targ or not sigma:
        return None
    # set data of first atom if all necessary data is present in the restraint
    ch2 = ch
    resn2 = resn
    altc2 = altc
    ins2 = ins
    atm2 = atm
    return ((ch1, resn1, ins1), altc1, atm1), \
           ((ch2, resn2, ins2), altc2, atm2), \
           (targ, sigma)

def parse_angle_restraint(restraint_line):
    """Parse the atoms, target and sigma from an angle restraint string.

    Return a tuple of (id_atom_1, id_atom_2, id_atom_3, (target, sigma),
    symm) where id_atom is a string tuple of (residue_id, alt_loc,
    atom_name). The middle atom is the vertex of the angle. symm is a tuple
    of three booleans, True for each atom that requests a symmetry operator
    ('symm Y').
    Return None if the angle restraint cannot be parsed.
    """
    if not restraint_line.startswith('exte angle'):
        return None
    is_next_chain, is_next_resnum, is_next_atom, is_next_ins = [False]*4
    is_next_altcode, is_next_symm, is_next_target = [False]*3
    is_next_sigma = False
    atoms = []
    symm = []
    ch, resn, atm, targ, sigma = [None]*5
    altc, ins = ' ', ' '
    atom_symm = False
    for word in restraint_line.split():
        if word == 'chain':
            is_next_chain = True
        elif word == 'resi' or word == 'residue':
            is_next_resnum = True
        elif word == 'atom':
            is_next_atom = True
        elif word == 'ins':
            is_next_ins = True
        elif word == 'alte':
            is_next_altcode = True
        elif word == 'symm':
            is_next_symm = True
        elif word == 'value':
            is_next_target = True
        elif word == 'sigma':
            is_next_sigma = True
        elif word == 'first' or word == 'next':
            # store the atom that was just read, if any
            if ch or resn is not None or atm:
                if not ch or resn is None or not atm:
                    return None
                atoms.append(((ch, resn, ins), altc, atm))
                symm.append(atom_symm)
            ch, resn, atm = None, None, None
            altc, ins = ' ', ' '
            atom_symm = False
        elif is_next_chain:
            is_next_chain = False
            ch = word
        elif is_next_resnum:
            is_next_resnum = False
            resn = int(word)
        elif is_next_atom:
            is_next_atom = False
            atm = word
        elif is_next_ins:
            is_next_ins = False
            if not word == '.':
                ins = word
        elif is_next_altcode:
            is_next_altcode = False
            if not word == '.':
                altc = word
        elif is_next_symm:
            is_next_symm = False
            if word.upper() == 'Y':
                atom_symm = True
        elif is_next_target:
            is_next_target = False
            targ = float(word)
        elif is_next_sigma:
            is_next_sigma = False
            sigma = float(word)
    # store the last atom
    if ch or resn is not None or atm:
        if not ch or resn is None or not atm:
            return None
        atoms.append(((ch, resn, ins), altc, atm))
        symm.append(atom_symm)
    # check if all three atoms have been set properly
    if len(atoms) != 3 or targ is None or sigma is None or sigma == 0:
        return None
    return atoms[0], atoms[1], atoms[2], (targ, sigma), tuple(symm)

def parse_atom_key(atom_words, column_dict):
    """Return the atom identity of pre-split mmCIF ATOM/HETATM words."""
    auth_asym_id = atom_words[column_dict['auth_asym_id']]
    auth_seq_id = int(atom_words[column_dict['auth_seq_id']])
    inscode = atom_words[column_dict['pdbx_PDB_ins_code']]
    if inscode == '?':
        inscode = ' '
    atom_name = atom_words[column_dict['auth_atom_id']]
    alt_code = atom_words[column_dict['label_alt_id']]
    if alt_code == '.':
        alt_code = ' '
    return ((auth_asym_id, auth_seq_id, inscode), alt_code, atom_name)


# Unused: index_atoms calls parse_atom_key directly.
# def parse_atom_line(line, column_dict):
#     """Return atom identity tuple from a mmCIF ATOM/HETATM line."""
#     return parse_atom_key(line.split(), column_dict)


def read_restraints(restraint_file):
    """Return distance and angle restraint tuples.

    The restraint tuples are atom strings and target/sd, followed for angle
    restraints by the per-atom symmetry flags.
    """
    distance_pairs = []
    angle_pairs = []
    lines = read(restraint_file)
    for line in lines:
        line = line.rstrip('\n')
        atom_pair = parse_distance_restraint(line)
        if atom_pair:
            distance_pairs.append(atom_pair)
            continue
        angle_restraint = parse_angle_restraint(line)
        if angle_restraint:
            atom1, atom2, atom3, target_sigma, symm = angle_restraint
            angle_pairs.append((atom1, atom2, atom3, target_sigma, symm))

    LOG.debug('%d distance restraints read.', len(distance_pairs))
    LOG.debug('%d angle restraints read.', len(angle_pairs))
    return distance_pairs, angle_pairs

if __name__ == '__main__':
    _main()
