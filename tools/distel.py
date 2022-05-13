#!/usr/bin/python3

"""
  Version 0.07 2022-02-02

  Calculate dihedral and distance restraint violation statistics.
  List rmsZ and outliers.
  Optionally create a YASARA macro to visualise restraint violations.

  Written by Robbie Joosten, Bart van Beusekom & Wouter Touw
  E-mail:    r.joosten@nki.nl, robie_joosten@hotmail.com

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


DISTANCE_MACRO = 'distel_distance_restraints.mcr'
DISTANCE_SCENE = 'distel_distance_restraints.sce'
TORSION_MACRO = 'distel_torsion_restraints.mcr'
TORSION_SCENE = 'distel_torsion_restraints.sce'
VERSION = 0.05

LOG = logging.getLogger(__name__)


def _main():
    """Compare restraints and structure.

    Restraint files must have the external restraints format (REFMAC).
    Distance and torsion angle restraints have been implemented.
    """
    fmt_version = 'distel (version {})'
    descr = '{} - {}'.format(fmt_version.format(VERSION),
                             'Analyze and visualize restraints.')
    parser = DefaultHelpParser(description=descr)
    parser.add_argument('-v', '--verbose', help='Verbose mode',
                        action='store_true')
    parser.add_argument('--version', help='Print version',
                        action=VersionActionStdOut,
                        version=fmt_version.format(VERSION))
    parser.add_argument('-m', '--yasmcr',
                        help='Write YASARA macro(s)? ({} and/or {})'.format(
                            DISTANCE_MACRO, TORSION_MACRO),
                        action='store_true')
    parser.set_defaults(verbose=False, print_version=False)
    parser.add_argument('pdb_file_path',
                        help='Path to input PDB file (.gz and .bz2 supported)',
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('rest_file_path',
                        help='Path to input restraint file '
                             '(.gz and .bz2 supported)',
                        type=lambda x: is_valid_file(parser, x,
                                                     empty_allowed=True))
    args = parser.parse_args()

    setup_logger(LOG, args.verbose)

    restraint_pairs, restraint_quads = read_restraints(args.rest_file_path)
    restraints = restraint_pairs + restraint_quads
    pdb_atoms = read_pdb(args.pdb_file_path)
    distance_matches, dihedral_matches = match_restraints_to_atoms(restraints,
                                                                   pdb_atoms)
    yasara_restraints = calc_distance_deviations(distance_matches, args.yasmcr)
    write_yasara_macro(yasara_restraints, DISTANCE_MACRO,
                       args.pdb_file_path, DISTANCE_SCENE)

    yasara_restraints = calc_dihedral_angle_deviations(dihedral_matches,
                                                       args.yasmcr)
    write_yasara_macro(yasara_restraints, TORSION_MACRO,
                       args.pdb_file_path, TORSION_SCENE, side_chains=False)

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
    

def calc_angle(v1, v2, v3):
    """Return angle in degrees."""
    u = v1 - v2
    w = v3 - v2
    return u.angle(w)*180/math.pi


def calc_angle_difference(target, found):
    """Return signed angle difference between target and found angle.

    Degrees
    """
    return (target - found + 180) % 360 - 180


def calc_angle_z(target, sigma, found):
    """Return Z in degrees."""
    return calc_angle_difference(target, found)/sigma


def calc_dihedral_angle(v1, v2, v3, v4):
    """Return dihedral angle in degrees."""
    ab = v1 - v2
    cb = v3 - v2
    db = v4 - v3
    # The dihedral angle is the angle between the planes spanned
    # by vectors ab/cb and db/cb
    u = ab ** cb
    v = db ** cb
    angle = u.angle(v)

    # Calculate the direction of normal vector of the plane spanned by
    # the two normal vectors above
    w = u ** v
    try:
        if cb.angle(w) > 0.001:
            angle = -angle
    except ZeroDivisionError:
        # The two normal vectors are parallel
        pass
    return angle*180/math.pi


def calc_dihedral_angle_deviations(matches, do_yasara):
    """Calculate rmsZ and outlier statistics for dihedral angle restraints.

    Return None if there are no matches or all restraints are ignored.
    """
    yasara_lines = []
    if len(matches) == 0:
        return None
    total_z_sq, num_restraints_ignored = 0.0, 0
    for match in matches:
        target = match[4][0]
        sigma = match[4][1]
        v1 = Vector.from_atom_record(match[0])
        v2 = Vector.from_atom_record(match[1])
        v3 = Vector.from_atom_record(match[2])
        v4 = Vector.from_atom_record(match[3])
        dih = calc_dihedral_angle(v1, v2, v3, v4)
        z = calc_angle_z(target, sigma, dih)
        total_z_sq += z**2
        msg = '{0:s} -- {1:s} -- {2:s} -- {3:s} | {4:7.2f} {5:7.2f} ' \
              '{6:7.2f} {7:7.2f}'.format(match[0][12:27], match[1][12:27],
                                         match[2][12:27], match[3][12:27],
                                         target, sigma, dih, z)
        if abs(z) > 4:
            LOG.info('OUTLIER: %s', msg)
        else:
            LOG.debug('         %s', msg)

        # generate line for YASARA plane
        if do_yasara:
            yasara_lines.append(yasara_restraint_plane(match,
                                                       abs(z)))
    if len(matches) == num_restraints_ignored:
        LOG.info('All restraints ignored')
        return None
    rmsz = math.sqrt(total_z_sq / (len(matches) - num_restraints_ignored))
    LOG.info('Torsion angle restraint rmsZ: %6.2f', rmsz)
    return yasara_lines


def calc_distance(v1, v2):
    """Return Euclidian distance between two vectors."""
    dif = v1 - v2
    return dif.l2_norm()


def calc_distance_deviations(matches, do_yasara):
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
        dist = calc_distance(Vector.from_atom_record(match[0]),
                             Vector.from_atom_record(match[1]))
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
                    Vector.from_atom_record(other_match[0]),
                    Vector.from_atom_record(other_match[1]))
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

            # generate line for YASARA arrow
            if do_yasara:
                yasara_lines.append(yasara_restraint_arrow(match,
                                                           violation))
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


def match_restraints_to_atoms(restraint_tuples, pdb_atoms):
    """Determine ATOM/HETATM lines that match restraint atom tuples.

    Return a tuple (distance, torsion) of lists of restraint target and sigma.
    """
    match1, match2, match3, match4 = None, None, None, None
    two_matches = []
    four_matches = []
    for restraint in restraint_tuples:
        for atom_line in pdb_atoms:
            if atom_line[21:27] == restraint[0][0] and \
                    atom_line[12:16].strip() == restraint[0][2] and \
                    atom_line[16] == restraint[0][1]:
                match1 = atom_line
            elif atom_line[21:27] == restraint[1][0] and \
                    atom_line[12:16].strip() == restraint[1][2] and \
                    atom_line[16] == restraint[1][1]:
                match2 = atom_line
            if len(restraint) == 5:
                if atom_line[21:27] == restraint[2][0] and \
                        atom_line[12:16].strip() == restraint[2][2] and \
                        atom_line[16] == restraint[2][1]:
                    match3 = atom_line
                elif atom_line[21:27] == restraint[3][0] and \
                        atom_line[12:16].strip() == restraint[3][2] and \
                        atom_line[16] == restraint[3][1]:
                    match4 = atom_line
        if match1 and match2 and (match3, match4) == (None, None):
            two_matches.append((match1, match2, restraint[2]))
        elif match1 and match2 and match3 and match4:
            four_matches.append((match1, match2, match3, match4, restraint[4]))
        match1, match2, match3, match4 = None, None, None, None

    LOG.debug('%d restraints matched with two atoms from PDB data',
              len(two_matches))
    LOG.debug('%d restraints matched with four atoms from PDB data',
              len(four_matches))
    return two_matches, four_matches


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
    return '{0:1s}{1:>4d}{2:1s}'.format(atom_words[1],
                                        int(atom_words[3]),
                                        ins), alte, atom_words[7]


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
    return ('{0:1s}{1:>4d}{2:1s}'.format(ch1, resn1, ins1), altc1, atm1), \
           ('{0:1s}{1:>4d}{2:1s}'.format(ch2, resn2, ins2), altc2, atm2), \
           (targ, sigma)


def parse_torsion_restraint(restraint_line):
    """Parse the atoms, target and sigma from torsion restraint string.

    Return a tuple of
    (id_atom_1, id_atom_2, id_atom_3, id_atom_4, (target, sigma)) where
    id_atom is a string tuple of (residue_id, alt_loc, atom_name).
    Return None if the torsion angle restraint cannot be parsed.

    A 'complete' external restraints format is expected
    'exte torsion first {0:s} next {1:s} next {2:s} next {3:s}'
    ' value {4:8.3f} sigma {5:8.3f} period 1'
    atoms are
    'chain {0:1s} resi {1:4d} inse {2:1s} atom {3:>4s} alte {4:1s}'

    The insertion codes and alternate codes are always expected.
    Raise a ValueError if insertion code, alternate code, value, sigma and
    period are not at the expected positions.
    Raise a ValueError if the period is not equal to 1.
    """
    restraint = None
    if not restraint_line.startswith('exte tors'):
        return None
    s = restraint_line.split()
    if len(s) != 52 or \
       (s[7], s[18], s[29], s[40]) != ('inse', )*4 or \
       (s[11], s[22], s[33], s[44]) != ('alte', )*4 or \
       (s[46], s[48], s[50]) != ('value', 'sigma', 'period'):
        raise ValueError('Unexpected restraint format: {}'.format(
            restraint_line))
    at1 = parse_complete_atom(s[3:13])
    at2 = parse_complete_atom(s[14:24])
    at3 = parse_complete_atom(s[25:35])
    at4 = parse_complete_atom(s[36:46])
    target, sigma = float(s[47]), float(s[49])
    period = int(s[51])
    if period != 1:
        raise ValueError('Period must be 1')
    restraint = (at1, at2, at3, at4, (target, sigma))
    return restraint


def read_pdb(pdb_file):
    """Return a list of ATOM and HETATM lines from PDB file."""
    lines = read(pdb_file)
    atoms = [l.strip() for l in lines if
             l.startswith('ATOM') or l.startswith('HETATM')]
    LOG.debug('%d atoms read.', len(atoms))
    return atoms


def read_restraints(restraint_file):
    """Return restraint tuple of distance and torsion restraint tuples.

    The restraint tuples are atom strings and target/sd.
    """
    distance_pairs = []
    torsion_quads = []
    lines = read(restraint_file)
    for line in lines:
        atom_pair = parse_distance_restraint(line.rstrip('\n'))
        if atom_pair:
            distance_pairs.append(atom_pair)
            continue
        atom_quad = parse_torsion_restraint(line.rstrip('\n'))
        if atom_quad:
            torsion_quads.append(atom_quad)
    LOG.debug('%d distance restraints read.', len(distance_pairs))
    LOG.debug('%d torsion angle restraints read.', len(torsion_quads))
    return distance_pairs, torsion_quads


def write_yasara_macro(yasara_restraints, yasara_file, pdb_file,
                       yasara_scene, side_chains=True):
    """Write YASARA macro for the generation of distance or torsion restraint
    scenes.
    """

    if not yasara_restraints or len(yasara_restraints) == 0:
        return None

    with open(yasara_file, 'w') as ym:
        ym.write('OnError Exit\n')
        ym.write('Fog 0\n')
        ym.write('ColorBG white\n')
        ym.write('Console off\n')
        ym.write('LoadPDB {}\n'.format(pdb_file))
        ym.write('ColorAtom all, gray\n')
        ym.write('Style Stick\n')
        if not side_chains:
            ym.write('HideAtom Protein Sidechain\n')
        for line in yasara_restraints:
            ym.write('{}\n'.format(line))
        ym.write('CenterAll\n')
        ym.write('NiceOriAll\n')
        ym.write('Zoom Steps=0\n')
        ym.write('SaveSce {}\n'.format(yasara_scene))
        ym.write('Exit\n')


def yasara_atom_selection(atom_record):
    """Return YASARA atom selection string from PDB format ATOM records."""
    atom_fmt = '{}{} Res {} Mol {}'
    alt = ''
    if atom_record[16] != ' ':
        alt = ' AltLoc {}'.format(atom_record[16])
    return atom_fmt.format(atom_record[12:16], alt, atom_record[22:27],
                           atom_record[21])


def yasara_restraint_arrow(match, violation):
    """Return ShowArrow YANACONDA command string."""
    y_fmt = 'ShowArrow Start=AtAtom, {}, End=AtAtom, {}, ' \
            'Radius=0.1, Color={}'
    return y_fmt.format(yasara_atom_selection(match[0]),
                        yasara_atom_selection(match[1]),
                        yasara_violation_colour(violation))


def yasara_restraint_plane(match, violation):
    """Return ShowPolygon YANACONDA command(s) string.

    Return Ci-1 - Ni - CAi for phi and also Ni - CAi - Ci
    Return CAi - Ci - Ni+1 for psi
    """
    col_fmt = 'ShowPolygon Coordinates=Atoms, Color={},  Alpha=70, ' \
              'Vertices=3, {}, {}, {}'
    share_fmt = 'ShowPolygon Coordinates=Atoms, Color=Gray, Alpha=50, ' \
                'Vertices=3, {}, {}, {}'
    if match[0][13] == 'C':
        # match C-N-CA-C
        phi_plane = col_fmt.format(yasara_violation_colour(violation),
                                   yasara_atom_selection(match[0]),
                                   yasara_atom_selection(match[1]),
                                   yasara_atom_selection(match[2]))
        share_plane = share_fmt.format(yasara_atom_selection(match[1]),
                                       yasara_atom_selection(match[2]),
                                       yasara_atom_selection(match[3]))
        return '{}\n{}'.format(phi_plane, share_plane)
    # match N-CA-C-N
    return col_fmt.format(yasara_violation_colour(violation),
                          yasara_atom_selection(match[1]),
                          yasara_atom_selection(match[2]),
                          yasara_atom_selection(match[3]))


def yasara_violation_colour(violation):
    """Map violation to integer between 120 and 359."""
    colour = int(359-24*violation)
    return colour if colour >= 120 else 120


if __name__ == '__main__':
    _main()
