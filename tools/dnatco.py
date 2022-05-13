#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests, sys, os, argparse, json, re, shutil, csv, gzip

# machine for POST request
# blackbox server, bigger machine, devel version (Mar 2020, identical to dnatco.org, both in v 3.2)
#dnatco_root='https://blackbox.ibt.biocev.org'
#dnatco_url=dnatco_root + '/devel'
#json2edits='/cgi-bin/json2edits_34devel_softer_chi.py'

# machine for POST request
# dnatco.org, still a bit small machine
dnatco_root='https://dnatco.datmos.org'
dnatco_url=dnatco_root + '/'
json2edits='/cgi-bin/json2edits_34_softer_chi.py'

# php file, most probably conformers_cif.php
php_file='conformers_cif.php'

# construct command line argument parser
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# define arguments
parser.add_argument("coordsfile", help="PDB or mmCIF filename")
parser.add_argument("-s", "--sigmafactor", default=1.0, help="Sigma scaling factor, smaller values produce stronger restraints.", type=float)
parser.add_argument("-c", "--cartesiancutoff", default=1.0, help="Cartesian RMSD cutoff in A, no restraints will be produced for steps with higher RMSD.", type=float)

args=parser.parse_args()

try:
  coordsfile = args.coordsfile
  sigmafactor = args.sigmafactor
  cartesiancutoff = args.cartesiancutoff

  print("Uploading and processing " + os.path.basename(coordsfile) +  " file ...")

  with open(coordsfile, 'rb') as f:
    r = requests.post(dnatco_url + '/' + php_file, files={'file': f}, data={'coordsfile': os.path.basename(coordsfile), 'sigmafactor': sigmafactor, 'cartesiancutoff': cartesiancutoff, 'output': 'json'})
   
  with open(os.path.basename(coordsfile) + ".html", 'wb') as f:
    f.write(r.text.encode('utf-8'))

  jsondata=json.loads(re.findall("^.*csvpath.*$",r.text,re.MULTILINE)[0])

 # with gzip.open(coordsfile + ".json.gz", 'wb') as jsoncsvfile:
 #   json.dump(jsondata,jsoncsvfile)

  print("Structure contains " + str(jsondata["jsonarr"]["num_std_residues"]) + " standard residues making " + str(jsondata["jsonarr"]["num_steps"]) + " steps.")
  print("Cartesian RMSD groups: " + str(jsondata["rmsd_counts"][0]) + "/" + str(jsondata["rmsd_counts"][1]) + "/" + str(jsondata["rmsd_counts"][2]))
  if (str(jsondata["jsonarr"]["num_models"]) != "1"):
    print("Found " + str(jsondata["jsonarr"]["num_models"]) + " MODELs in the structure, restraints will be created for MODEL 1 only." )

  #print("Preparing restraints edits for dinucleotide steps within " + str(cartesiancutoff) + "Å scaled by " + str(sigmafactor) + " factor ...")

  #with requests.post(dnatco_root + json2edits, json=jsondata["edits_dict"] ) as r:
    #edits_json=json.loads(r.text)
  #with open(os.path.basename(coordsfile) + "_s" + str(sigmafactor) + "_c" + str(cartesiancutoff) + "_edits", 'wb') as f:
    #for line in edits_json["edits_lines"]:
      #f.write(line + "\n")

# editspath so far points to unfiltered edits file, better version is above
#  for path in (jsondata["csvpath"], jsondata["editspath"], jsondata["jsonpath"]):
# restraints for Phenix, REFMAC, and MMB are work in progress
  #for path in (jsondata["csvpath"], jsondata["jsonpath"]):
    #filename = path.split('/')[-1]
    #print("Downloading "  + filename  + " ...")
  with requests.get(dnatco_url + '/' + jsondata["jsonpath"], stream=True) as r:
    with open(os.path.basename(coordsfile) + ".dnatco.json.gz", 'wb') as f:
      shutil.copyfileobj(r.raw, f)

  #jsoncsv = {}
  #with open(jsondata["csvpath"].split('/')[-1]) as csvfile:
    #reader = csv.DictReader(csvfile)
    #for row in reader:
      #jsoncsv[row["step_ID"]] = row

  #with open(jsondata["csvpath"].split('/')[-1] + ".json", 'wb') as jsoncsvfile:
    #json.dump(jsoncsv,jsoncsvfile)

  print(os.path.basename(coordsfile) +  " done")

except:
  print("something went wrong")
  raise
