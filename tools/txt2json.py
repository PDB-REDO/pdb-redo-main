#!/usr/bin/python3

"""
Convert a data.txt file (first argument) to data.json.

o This script is compatible python3. 
o Only save the last data.txt line (without warning) if the input consists of
  multiple data lines from the same data.txt file or multiple concatenated
  data.txt files.
o Clean data.txt file:
    - Add NAs to trailing columns.
o Print a warning when the number of data.txt columns is
  larger than the expected number of columns defined in the dictionary
  JSON_VAR.
o Translate 'NA', '<NA>' and 'none' always to None (Python) -> null (JSON).
o JSON does not support date types. 'YYYY-mm-dd' TIME is a string here.
o Write nested JSON output with keys: 'pdbid' and 'properties'
"""

from __future__ import print_function

import argparse
from ast import Or
import json
import re
import sys

from collections import OrderedDict


__version__ = '2.1'

# In data.txt files the columns are always in the same order
# and don't have headers
# The position in JSON_VAR indicates the column type
JSON_VAR = OrderedDict()
JSON_VAR['PDBID'] = 'str'
JSON_VAR['VERSION'] = 'float'
JSON_VAR['RFACT'] = 'float'
JSON_VAR['RFREE'] = 'float'
JSON_VAR['RCAL'] = 'float'
JSON_VAR['RFCAL'] = 'float'
JSON_VAR['SIGRFCAL'] = 'float'
JSON_VAR['RFCALUNB'] = 'float'
JSON_VAR['RFCALZ'] = 'float'
JSON_VAR['RTLS'] = 'float'
JSON_VAR['RFTLS'] = 'float'
JSON_VAR['SIGRFTLS'] = 'float'
JSON_VAR['RFTLSUNB'] = 'float'
JSON_VAR['RFTLSZ'] = 'float'
JSON_VAR['RFIN'] = 'float'
JSON_VAR['RFFIN'] = 'float'
JSON_VAR['SIGRFFIN'] = 'float'
JSON_VAR['RFFINUNB'] = 'float'
JSON_VAR['RFFINZ'] = 'float'
JSON_VAR['RFRRAT'] = 'float'
JSON_VAR['BBEST'] = 'str'
JSON_VAR['TLSBEST'] = 'str'
JSON_VAR['BLTBEST'] = 'str'
JSON_VAR['NWATDEL'] = 'int'
JSON_VAR['NBBFLIP'] = 'int'
JSON_VAR['NSCBLT'] = 'int'
JSON_VAR['NDROTA'] = 'int'
JSON_VAR['HBFLIP'] = 'int'
JSON_VAR['STFLIP'] = 'int'
JSON_VAR['NCHIRFX'] = 'int'
JSON_VAR['OZPAK1'] = 'float'
JSON_VAR['NZPAK1'] = 'float'
JSON_VAR['FZPAK1'] = 'float'
JSON_VAR['OZPAK2'] = 'float'
JSON_VAR['NZPAK2'] = 'float'
JSON_VAR['FZPAK2'] = 'float'
JSON_VAR['OZRAMA'] = 'float'
JSON_VAR['NZRAMA'] = 'float'
JSON_VAR['FZRAMA'] = 'float'
JSON_VAR['OCHI12'] = 'float'
JSON_VAR['NCHI12'] = 'float'
JSON_VAR['FCHI12'] = 'float'
JSON_VAR['OBCONF'] = 'float'
JSON_VAR['NBCONF'] = 'float'
JSON_VAR['FBCONF'] = 'float'
JSON_VAR['OBRMSZ'] = 'float'
JSON_VAR['NBRMSZ'] = 'float'
JSON_VAR['FBRMSZ'] = 'float'
JSON_VAR['OARMSZ'] = 'float'
JSON_VAR['NARMSZ'] = 'float'
JSON_VAR['FARMSZ'] = 'float'
JSON_VAR['OBUMPS'] = 'int'
JSON_VAR['NBUMPS'] = 'int'
JSON_VAR['FBUMPS'] = 'int'
JSON_VAR['OHBUNS'] = 'int'
JSON_VAR['NHBUNS'] = 'int'
JSON_VAR['FHBUNS'] = 'int'
JSON_VAR['OGFOLD'] = 'float'
JSON_VAR['NGFOLD'] = 'float'
JSON_VAR['FGFOLD'] = 'float'
JSON_VAR['PROG'] = 'str'
JSON_VAR['DYEAR'] = 'int'
JSON_VAR['RESOLUTION'] = 'float'
JSON_VAR['DATARESH'] = 'float'
JSON_VAR['DATARESL'] = 'float'
JSON_VAR['NREFCNT'] = 'int'
JSON_VAR['TSTCNT'] = 'int'
JSON_VAR['TSTPRC'] = 'float'
JSON_VAR['NTSTCNT'] = 'int'
JSON_VAR['REFPATM'] = 'float'
JSON_VAR['AAXIS'] = 'float'
JSON_VAR['BAXIS'] = 'float'
JSON_VAR['CAXIS'] = 'float'
JSON_VAR['ALPHA'] = 'float'
JSON_VAR['BETA'] = 'float'
JSON_VAR['GAMMA'] = 'float'
JSON_VAR['BAVER'] = 'float'
JSON_VAR['BWILS'] = 'float'
JSON_VAR['BREFTYPE'] = 'str'
JSON_VAR['SOLVENT'] = 'str'
JSON_VAR['VDWPROBE'] = 'float'
JSON_VAR['IONPROBE'] = 'float'
JSON_VAR['RSHRINK'] = 'float'
JSON_VAR['DOTLS'] = 'bool'
JSON_VAR['NTLS'] = 'int'
JSON_VAR['OPTTLSG'] = 'str'
JSON_VAR['ORITLS'] = 'bool'
JSON_VAR['LEGACY'] = 'bool'
JSON_VAR['SPACEGROUP'] = 'str'
JSON_VAR['RSCCB'] = 'int'
JSON_VAR['RSCCW'] = 'int'
JSON_VAR['RSRB'] = 'int'
JSON_VAR['RSRW'] = 'int'
JSON_VAR['OWBMPS'] = 'float'
JSON_VAR['NWBMPS'] = 'float'
JSON_VAR['FWBMPS'] = 'float'
JSON_VAR['OHBSAT'] = 'float'
JSON_VAR['NHBSAT'] = 'float'
JSON_VAR['FHBSAT'] = 'float'
JSON_VAR['URESO'] = 'float'
JSON_VAR['CCWOLD'] = 'float'
JSON_VAR['CCWFIN'] = 'float'
JSON_VAR['ZCCW'] = 'float'
JSON_VAR['CCFOLD'] = 'float'
JSON_VAR['CCFFIN'] = 'float'
JSON_VAR['ZCCF'] = 'float'
JSON_VAR['WAVELENGTH'] = 'float'
JSON_VAR['ISTWIN'] = 'bool'
JSON_VAR['SOLVD'] = 'float'
JSON_VAR['EXPTYP'] = 'str'
JSON_VAR['COMPLETED'] = 'float'
JSON_VAR['NOPDB'] = 'bool'
JSON_VAR['NOSF'] = 'bool'
JSON_VAR['USIGMA'] = 'bool'
JSON_VAR['ZCALERR'] = 'bool'
JSON_VAR['TIME'] = 'YYYY-mm-dd'
JSON_VAR['RESOTYPE'] = 'int'
JSON_VAR['FALSETWIN'] = 'bool'
JSON_VAR['TOZRAMA'] = 'int'
JSON_VAR['TFZRAMA'] = 'int'
JSON_VAR['TOCHI12'] = 'int'
JSON_VAR['TFCHI12'] = 'int'
JSON_VAR['TOZPAK2'] = 'int'
JSON_VAR['TFZPAK2'] = 'int'
JSON_VAR['TOWBMPS'] = 'int'
JSON_VAR['TFWBMPS'] = 'int'
JSON_VAR['TOHBSAT'] = 'int'
JSON_VAR['TFHBSAT'] = 'int'
JSON_VAR['OHRMSZ'] = 'float'
JSON_VAR['NHRMSZ'] = 'float'
JSON_VAR['FHRMSZ'] = 'float'
JSON_VAR['TOZPAK1'] = 'int'
JSON_VAR['TFZPAK1'] = 'int'
JSON_VAR['NLOOPS'] = 'int'
JSON_VAR['NMETALREST2'] = 'int'
JSON_VAR['OSZRAMA'] = 'float'
JSON_VAR['NSZRAMA'] = 'float'
JSON_VAR['FSZRAMA'] = 'float'
JSON_VAR['OSCHI12'] = 'float'
JSON_VAR['NSCHI12'] = 'float'
JSON_VAR['FSCHI12'] = 'float'
JSON_VAR['SRFRRAT'] = 'float'
JSON_VAR['ZRFRRATCAL'] = 'float'
JSON_VAR['ZRFRRATTLS'] = 'float'
JSON_VAR['ZRFRRATFIN'] = 'float'
JSON_VAR['GOT_PROT'] = 'bool'
JSON_VAR['GOT_NUC'] = 'bool'
JSON_VAR['NNUCLEICREST'] = 'int'
JSON_VAR['OBPHBRMSZ'] = 'float'
JSON_VAR['NBPHBRMSZ'] = 'float'
JSON_VAR['FBPHBRMSZ'] = 'float'
JSON_VAR['OSHEAR'] = 'float'
JSON_VAR['NSHEAR'] = 'float'
JSON_VAR['FSHEAR'] = 'float'
JSON_VAR['OSTRETCH'] = 'float'
JSON_VAR['NSTRETCH'] = 'float'
JSON_VAR['FSTRETCH'] = 'float'
JSON_VAR['OBUCKLE'] = 'float'
JSON_VAR['NBUCKLE'] = 'float'
JSON_VAR['FBUCKLE'] = 'float'
JSON_VAR['OPROPEL'] = 'float'
JSON_VAR['NPROPEL'] = 'float'
JSON_VAR['FPROPEL'] = 'float'
JSON_VAR['OCONFAL'] = 'float'
JSON_VAR['NCONFAL'] = 'float'
JSON_VAR['FCONFAL'] = 'float'
JSON_VAR['TOCONFAL'] = 'int'
JSON_VAR['TNCONFAL'] = 'int'
JSON_VAR['TFCONFAL'] = 'int'
JSON_VAR['ODNRMSD'] = 'float'
JSON_VAR['NDNRMSD'] = 'float'
JSON_VAR['FDNRMSD'] = 'float'
JSON_VAR['OBPGRMSZ'] = 'float'
JSON_VAR['NBPGRMSZ'] = 'float'
JSON_VAR['FBPGRMSZ'] = 'float'
JSON_VAR['TOBPGRMSZ'] = 'int'
JSON_VAR['TFBPGRMSZ'] = 'int'
JSON_VAR['FSCWCAL'] = 'float'
JSON_VAR['FSCFCAL'] = 'float'
JSON_VAR['FSCWTLS'] = 'float'
JSON_VAR['FSCFTLS'] = 'float'
JSON_VAR['FSCWFIN'] = 'float'
JSON_VAR['FSCFFIN'] = 'float'

# The data is space-delimited, except for the spacegroup between quotes
RE_COL = re.compile(r"'[^']+'|[^'\s]+")


def add_na(data_line):
    """Add NA to trailing columns.

    Return a list of string columns.
    """
    # Extract data values
    match = RE_COL.findall(data_line)
    match_len = len(match)
    n_var = len(JSON_VAR)
    if match_len > n_var:
        # Assume first column is always PDBID
        pdb_id = match[0]
        warning(f"Data line for {pdb_id} longer ({match_len}) than expected ({n_var})")
    if match_len < n_var:
        match.extend(['NA']*(n_var - match_len))
    return match


def clean_data(fh_data_txt):
    """Clean data.txt.

    Currently only add NAs to trailing columns.

    Return a list of values for the data.txt line (or the last line in case
    multiple lines were given as input).
    """
    with fh_data_txt as handle:
        # Select data lines
        lines = [add_na(l.strip()) for l in handle if not l.startswith('#')]
    if len(lines) == 0:
        raise ValueError('No data lines found')
    # Return last line only
    return lines[-1]


def parse(clean_list):
    """Convert clean_list string list of values to python data types.

    Return a dict.
    """
    def cast(value, var_type):
        """Cast value to type var_type.

        Raise ValueError for unsupported formats.
        """
        if value in ['NA', 'none', '<NA>', 'NaN', '-nan']:
            return None
        if value[0] == "'" and value[-1] == "'":
            # Remove single quotes around spacegroup
            value = value[1:-1]
        if var_type in ['str', 'YYYY-mm-dd']:
            return value
        if var_type == 'float':
            return float(value)
        if var_type == 'int':
            return int(value)
        if var_type == 'bool':
            if value == 'T':
                return True
            elif value == 'F':
                return False
            else: 
                return bool(int(value))
        raise ValueError(f'Invalid type: {var_type}')

    return {var: cast(value, var_type) for (var, var_type), value in
            zip(JSON_VAR.items(), clean_list)}



def warning(*objs):
    """Print a warning to stderr."""
    print("WARNING: ", *objs, file=sys.stderr)


def write_json(parsed_data, fh_out):
    """Write the JSON file."""    
    nested_outp={}
    nested_outp['pdbid'] = parsed_data['PDBID']      #define pdbid as key in output json dict
    nested_outp['properties'] = parsed_data          #define properties as key in output json dict
    del nested_outp['properties']['PDBID']           #delete the PDB-ID key from the 'properties'
    
    with fh_out as handle:    
        json.dump(nested_outp, handle, sort_keys=True, indent=4)


def main():
    """Convert a PDB_REDO data.txt file to JSON."""
    parser = argparse.ArgumentParser(
        description=f'txt2json version {__version__}')
    parser.add_argument('-c', '--comment',
                        help='optional comment')
    parser.add_argument('-i', '--infile', nargs='?',
                        help='data.txt file location (defaults to stdin)',
                        type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', nargs='?',
                        help='output file location (defaults to stdout)',
                        type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args()

    # Show help if no input file has been given
    # via neither a pipe nor the infile argument
    if len(sys.argv) < 2:
        parser.parse_args(['-h'])
        sys.exit(1)

    data = clean_data(args.infile)
    parsed_data = parse(data)
    parsed_data['COMMENT'] = args.comment
    write_json(parsed_data, args.outfile)


if __name__ == '__main__':
    main()
