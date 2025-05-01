#!/usr/bin/python

"""
  Version 0.08 2025-05-01

  Calculate  distance restraint violation statistics.
  List rmsZ and outliers.
  Optionally create a YASARA macro to visualise restraint violations.

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

import logging
import math
import numpy as np

from helpers import DefaultHelpParser, Vector, VersionActionStdOut, \
                    is_valid_file, read, setup_logger

VERSION = 0.08

LOG = logging.getLogger(__name__)


def _main():
    """Compare restraints and structure.

    Restraint files must have the external restraints format (REFMAC).
    Distance and torsion angle restraints have been implemented.
    """
    fmt_version = 'distel (version {})'
    descr = '{} - {}'.format(fmt_version.format(VERSION),
                             'Analyze restraints.')
    parser = DefaultHelpParser(description=descr)
    parser.add_argument('-v', '--verbose', help='Verbose mode',
                        action='store_true')
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

    restraint_pairs = read_restraints(args.rest_file_path)
    restraints = restraint_pairs 
    mmcif_atoms, atom_site_columns = read_mmcif(args.mmcif_file_path)
    distance_matches = match_restraints_to_atoms(restraints,mmcif_atoms,atom_site_columns)

    calc_distance_deviations(distance_matches,atom_site_columns)
    
    
def jackknife(x, func):
    """Jackknife estimate of the estimator func"""
    n = len(x)
    idx = np.arange(n)
    return sum(func(x[idx!=i]) for i in range(n))/float(n)

def jackknife_stdev(x, func):
    """Jackknife estimate of the variance of the estimator func."""
    n = len(x)
    idx = np.arange(n)
    j_est = jackknife(x, func)
    return np.sqrt((n-1)/(n + 0.0) * sum((func(x[idx!=i]) - j_est)**2.0 for i in range(n)))

def rms(x):
    return np.sqrt(np.mean(x**2))
    

def calc_distance(v1, v2):
    """Return Euclidian distance between two vectors."""
    dif = v1 - v2
    return dif.l2_norm()


def calc_distance_deviations(matches, atom_site_columns):
    """Calculate rmsZ and outlier statistics for distance restraints.

    Return lines for YASARA macro if do_yasara is True.
    Return None if there are no matches or all restraints are ignored.
    """
    if len(matches) == 0:
        LOG.info('Distance restraint rmsZ: NA')
        return None
    total_z_sq, num_restraints_ignored = 0.0, 0
    yasara_lines = []
    lzlength = []
    for match in matches:
        target = match[2][0]
        sigma = match[2][1]
        dist = calc_distance(Vector.from_atom_record(match[0],atom_site_columns),
                             Vector.from_atom_record(match[1],atom_site_columns))
        z = (dist - target) / sigma
        violation = abs(z)
        # check if there is not another restraint for the same atoms with a
        # target closer to reality
        keep_match = True
        for other_match in matches:
            if match[0] == other_match[0] and match[1] == other_match[1]:
                # if getting here, there are indeed two restraints for the same
                # atoms. Filter out the target with greatest deviation
                other_target = other_match[2][0]
                other_sigma = other_match[2][1]
                other_dist = calc_distance(
                    Vector.from_atom_record(other_match[0], atom_site_columns),
                    Vector.from_atom_record(other_match[1],atom_site_columns))
                other_z = (other_dist - other_target) / other_sigma
                other_violation = abs(other_z)
                if other_violation < violation:
                    keep_match = False
        if keep_match:
            lzlength.append(z)
            total_z_sq += z**2
            msg = '{0:s} -- {1:s} | {2:6.2f} {3:6.2f} {4:6.2f} ' \
                  '{5:6.2f}'.format(match[0][12:27], match[1][12:27], target,
                                    sigma, dist, z)
            if abs(z) > 4:
                LOG.info('OUTLIER: %s', msg)
            else:
                LOG.debug('         %s', msg)
        else:
            num_restraints_ignored += 1
            
    azlength = np.array(lzlength)
    #jrmsz = jackknife(azlength, rms)
    stdev = jackknife_stdev(azlength, rms)
    #print (jrmsz, stdev)
    if len(matches) == num_restraints_ignored:
        LOG.info('All restraints ignored')
        LOG.info('Distance restraint rmsZ: NA')
        return None
    rmsz = math.sqrt(total_z_sq / (len(matches) - num_restraints_ignored))
    LOG.info('Distance restraint rmsZ: %6.3f %6.3f', rmsz, stdev)
    return yasara_lines


def match_restraints_to_atoms(restraint_tuples, mmcif_atoms, atom_site_columns):
    """Determine ATOM/HETATM lines that match restraint atom tuples.

    Return a tuple (distance, torsion) of lists of restraint target and sigma.
    """
    
    column_dict = {col:i for i, col in enumerate(atom_site_columns)}
    
   
    
    match1, match2 = None, None
    two_matches = []
    for restraint in restraint_tuples:
        for atom_line in mmcif_atoms:
            
            atom_info = parse_atom_line(atom_line,column_dict)
            #print(atom_info, restraint)
            if atom_info == restraint[0]:
                match1 = atom_line
            elif atom_info == restraint[1]:
                match2 = atom_line
        if match1 and match2:
            two_matches.append((match1, match2, restraint[2]))
        match1, match2 = None, None

    LOG.debug('%d restraints matched with two atoms from model data',
              len(two_matches))
    return two_matches

def read_mmcif(mmcif_file):
    """Return a list of ATOM and HETATM lines from model file."""
    lines = read(mmcif_file)
    atom_site_columns = []
    atoms = []
    capture = False
    previous_line_was_loop = False
    
    for line in lines:
        line = line.strip()
        if previous_line_was_loop and line.startswith('_atom_site'):
            capture = True
            atom_site_columns.append(line.split('.')[1])  # Add the current line
        elif capture:
            if line.startswith('_atom_site'):
                atom_site_columns.append(line.split('.')[1])  # Continue adding
            else:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atoms.append(line)
        
        # Check if current line is 'loop_' for the next iteration
        previous_line_was_loop = line == 'loop_'
        
    
    LOG.debug('%d atoms read.', len(atoms))
    return atoms, atom_site_columns

def parse_complete_atom(atom_words):
    """Parse the list atom_words in external restraint format.

    Return a tuple (residue_id, alt_loc, atom_name)
    Return a ValueError if the list is not of the expected length or if the
    a string cannot be cast to the expected type.
    """
    exp = 'chain <chain> resi <resi> inse <inse> atom <atom> alte <alte>'
    if not len(atom_words) == 10:
        raise ValueError('Expected atom format: {}'.format(exp))
    ins, alte = atom_words[5], atom_words[9]
    ins = ' ' if ins == '.' else ins
    alte = ' ' if alte == '.' else alte
    return (atom_words[1],int(atom_words[3]),ins), alte, atom_words[7]


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

def parse_atom_line(line, column_dict):
    
    auth_asym_id = line.split()[column_dict['auth_asym_id']]
    auth_seq_id = int(line.split()[column_dict['auth_seq_id']])
    inscode = ' ' if line.split()[column_dict['pdbx_PDB_ins_code']] == '?' else line.split()[column_dict['pdbx_PDB_ins_code']]
    atom_name = line.split()[column_dict['auth_atom_id']]
    alt_code =  ' ' if line.split()[4] == '.' else line.split()[4]
     
    atom_info = ((auth_asym_id, auth_seq_id, inscode), alt_code, atom_name)
    
    return atom_info


def read_restraints(restraint_file):
    """Return restraint tuple of distance and torsion restraint tuples.

    The restraint tuples are atom strings and target/sd.
    """
    distance_pairs = []
    lines = read(restraint_file)
    for line in lines:
        atom_pair = parse_distance_restraint(line.rstrip('\n'))
        if atom_pair:
            distance_pairs.append(atom_pair)
            continue

    LOG.debug('%d distance restraints read.', len(distance_pairs))
    return distance_pairs

if __name__ == '__main__':
    _main()
