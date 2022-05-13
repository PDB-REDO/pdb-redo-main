import os
import re
import sys
import json


def dssr_data(dssr_file):
    with open(dssr_file) as file:
        data = json.load(file)
        DNA = ["DA", "DC", "DG", "DT", "DU"]
        DNAnc = ["08Q","0AD","0DC","0DG","0DT","1AP","1CC","1W5","1WA","2AR","2AT","2BD","2BT","2BU","2DA","2EG","2GT","2NT","2OT","2PR","2ST","3DR","47C","4PC","4PD","4PE","5CM","5FC","5HC","5HU","5IU","5NC","5PC","5PY","5SE","64T","6HA","6HC","6HG","6HT","6MA","6OG","7DA","7GU","85Y","8FG","8MG","8OG","A5L","ABR","ABS","AF2","AFG","ASU","B7C","BGM","BOE","BRU","C34","C6G","C7R","C7S","CAR","CB2","CBR","CDW","CFZ","CTG","D3","D33","DA","DC","DDG","DDN","DFT","DG","DI","DN","DNR","DOC","DRZ","DT","DU","DUZ","DXD","DZM","E","E1X","EDA","EDC","EIT","FDG","FMG","G35","G36","GF2","GMU","GN7","GNE","GSR","GSS","HEU","HOL","IGU","IMC","LCG","LGP","M1G","MA7","MBZ","MG1","MRG","MTR","NCU","NCX","NDU","NMS","NMT","NTT","P","P2T","P2U","PDU","PGN","PPW","PRN","R","S2M","S6G","SOS","T2T","T49","T4S","T5O","T5S","TA3","TAF","TC1","TDY","TFE","TLC","TLN","TTM","U2N","UFT","UMS","US3","USM","X","X4A","XAE","XAL","XCL","XCR","XCS","XCY","XGA","XGL","XGR","XTL","XTR","XTY","XUA","XUG","Y","YRR","Z","ZDU"]
        RNA = ["A", "C", "G", "T", "U"]
        RNAnc = ["0A","0C","0G","0U","0U1","125","126","127","12A","1MA","1MG","23G","2AU","2MA","2MG","2MU","2OM","4AC","4OC","4SU","5BU","5CF","5IC","5MC","5MU","6IA","6NW","6OP","70U","7MG","A","A23","A2L","A2M","A44","A5M","A5O","A6A","A6C","A6G","A6U","AET","AP7","AVC","BMP","C","C2L","C43","CCC","CH","F2T","FHU","G","G25","G2L","G48","G7M","GAO","H2U","I","IC","IG","IU","KAG","LCC","M2G","M5M","MA6","MGV","MIA","MNU","N","N5M","N6G","NF2","O2G","OMC","OMG","OMU","ONE","P5P","PGP","PQ1","PSU","PYO","QUO","RSQ","S4C","SMT","SSU","SUR","T","T6A","U","U2L","U36","U8U","UAR","UOB","UR3","YG","YYG"] 
        try: 
            for p in data["pairs"]:
                type_bond = p["name"]
                if type_bond == "WC" or type_bond == "Wobble":
                    nucleotides = p["bp"]
                    chain1 = p["nt1"][0:1]
                    try:
                        pos = p["nt1"].index("^")   
                        ins1 = p["nt1"][pos+1:pos+2]
                    except:
                        ins1 = "."
                    m = re.search("(-?\d+)", p["nt1"][2:])
                    res1 = m.string[:m.start()]
                    if res1 in DNA:
                        typ1 = "DNA"
                    elif res1 in RNA:
                        typ1 = "RNA"
                    elif res1 in DNAnc:
                        typ1 = "DNA"
                    elif res1 in RNAnc:
                        typ1 = "RNA"                        
                    else:
                        typ1 = "other"   
                    seqnum1 = re.findall("(-?\d+)", p["nt1"][2:])[-1]                
                    chain2 = p["nt2"][0:1]
                    try:
                        pos = p["nt2"].index("^")
                        ins2 = p["nt2"][pos+1:pos+2]
                    except:  
                        ins2 = "."
                    m = re.search("(-?\d+)", p["nt2"][2:])
                    res2 = m.string[:m.start()]
                    if res2 in DNA:
                        typ2 = "DNA"
                    elif res2 in RNA:
                        typ2 = "RNA"
                    elif res2 in DNAnc:
                        typ2 = "DNA"
                    elif res2 in RNAnc:
                        typ2 = "DNA" 
                    else:
                        typ2 = "other"  
                    seqnum2 = re.findall("(-?\d+)", p["nt2"][2:])[-1]    
                    a =  p["hbonds_desc"].split(",")
                    for i in a:
                        if "O6" in i and "N4" in i:
                            if i[0:2] == "O6":
                                bond = "O6-N4"
                            else:
                                bond = "N4-O6"
                        elif "N1" in i and "N3" in i:
                            if i[0:2] == "N1":
                                bond = "N1-N3"
                            else:
                                bond = "N3-N1"
                        elif "O2" in i and "N2" in i:
                            if i[0:2] == "O2":
                                bond = "O2-N2"
                            else:
                                bond = "N2-O2"
                        elif "O4" in i and "N6" in i:
                            if i[0:2] == "O4":
                                bond = "O4-N6"
                            else:
                                bond = "N6-O4"
                        elif "O2" in i and "N1" in i:
                            if i[0:2] == "O2":
                                bond = "O2-N1"
                            else:
                                bond = "N1-O2"
                        elif "O6" in i and "N3" in i:
                            if i[0:2] == "O6":
                                bond = "O6-N3"
                            else:
                                bond = "N3-O6"
                        else:
                            bond = "other"

                        if type_bond == "WC" and bond != "other":
                            if "O6" in bond and "N4" in bond:
                                if "c" in nucleotides.lower() and "g" in nucleotides.lower():
                                    if typ1 == "DNA" and typ2 == "DNA":
                                        bond_mean = 2.901
                                        bond_std = 0.095
                                    elif typ1 == "RNA" and typ2 == "RNA":
                                        bond_mean = 2.916
                                        bond_std = 0.088
                                    else:
                                        bond_mean = 2.889
                                        bond_std = 0.051
                            elif "N1" in bond and "N3" in bond:
                                if "a" in nucleotides.lower() and "t" in nucleotides.lower():
                                    if typ1 == "DNA" and typ2 == "DNA":
                                        bond_mean = 2.825
                                        bond_std = 0.053
                                    else:  
                                    #RNA-RNA and RNA-DNA pairs are treated the same as RNA-RNA pairs are very rare
                                        bond_mean = 2.800
                                        bond_std = 0.029
                                elif "a" in nucleotides.lower() and "u" in nucleotides.lower():
                                    if typ1 == "RNA" and typ2 == "RNA":
                                        bond_mean = 2.829
                                        bond_std = 0.051
                                    else:
                                    #DNA-DNA and RNA-DNA pairs are treated the same as DNA-DNA pairs are very rare
                                        bond_mean = 2.799
                                        bond_std = 0.025
                                elif "c" in nucleotides.lower() and "g" in nucleotides.lower():
                                    if typ1 == "DNA" and typ2 == "DNA":
                                        bond_mean = 2.907
                                        bond_std = 0.055
                                    elif typ1 == "RNA" and typ2 == "RNA":
                                        bond_mean = 2.901
                                        bond_std = 0.049
                                    else:
                                        bond_mean = 2.875
                                        bond_std = 0.053
                            elif "O2" in bond and "N2" in bond:
                                if "c" in nucleotides.lower() and "g" in nucleotides.lower():
                                    if typ1 == "DNA" and typ2 == "DNA":
                                        bond_mean = 2.830
                                        bond_std = 0.078
                                    elif typ1 == "RNA" and typ2 == "RNA":
                                        bond_mean = 2.819
                                        bond_std = 0.070
                                    else:
                                        bond_mean = 2.779
                                        bond_std = 0.098
                            elif "O4" in bond and "N6" in bond:
                                if "a" in nucleotides.lower() and "t" in nucleotides.lower():
                                    if typ1 == "DNA" and typ2 == "DNA":
                                        bond_mean = 2.999
                                        bond_std = 0.099
                                    else:  
                                    #RNA-RNA and RNA-DNA pairs are treated the same as RNA-RNA pairs are very rare
                                        bond_mean = 2.939
                                        bond_std = 0.058
                                elif "a" in nucleotides.lower() and "u" in nucleotides.lower():    
                                    if typ1 == "RNA" and typ2 == "RNA":
                                        bond_mean = 2.974
                                        bond_std = 0.094
                                    else:
                                    #DNA-DNA and RNA-DNA pairs are treated the same as DNA-DNA pairs are very rare
                                        bond_mean = 2.963
                                        bond_std = 0.039
                            else:
                                continue
                                
                            try:
                                string = "exte dist first chain " + chain1 + " resi " + seqnum1 + " ins " + ins1 + " atom " + bond[:2] + " second chain " + chain2 + " resi " + seqnum2 + " ins " + ins2 + " atom " + bond[3:5] + " value " + str(bond_mean) + " sigma " + str(bond_std) + " type 1"
                                print(string)
                            except:
                                UnboundLocalError
                        #elif type_bond == "Wobble" and bond != "other":
                            #if "O2" in bond and "N1" in bond:
                                #if "u" in nucleotides.lower() and "g" in nucleotides.lower():
                                    #bond_mean = 2.820
                                    #bond_std = 0.
                            #elif "O6" in bond and "N3" in bond:
                                #if "u" in nucleotides.lower() and "g" in nucleotides.lower():
                                    #bond_mean = 2.826
                                    #bond_std = 0.109
                            #try:
                                #string = "exte dist first chain " + chain1 + " resi " + seqnum1 + " ins " + ins1 + " atom " + bond[:2] + " second chain " + chain2 + " resi " + seqnum2 + " ins " + ins2 + " atom " + bond[3:5] + " value " + str(bond_mean) + " sigma " + str(bond_std) + " type 1"
                                #print(string)
                            #except:
                                #UnboundLocalError
        except: 
            NameError
        file.close()
    return()
        
###INS ALTIJD .
###TYPE 1
            
def run_dssr(pdb_file, dssr, output_dir):
    os.mkdir(output_dir)
    path = pdb_file.split("/")
    file_name = path[-1]    
    out_name = (file_name + ".json")
    out_file = output_dir + "/" + out_name
#    print(out_file)
    search_command = (dssr + " -i=" + pdb_file +  " -o=" + out_file + " --json"  + " --more" + " --non-pair > /dev/null 2>&1 ")  
    os.system(search_command)
    return(out_file)    

pdb_file = os.path.abspath(sys.argv[1])
dssr_path = os.path.abspath(sys.argv[2])
output_path = "dssr_temp"
dssr_file = run_dssr(pdb_file, dssr_path, output_path)
dssr_data(dssr_file)

#Clean up
search_command = ('rm -r dssr_temp')
os.system(search_command)

