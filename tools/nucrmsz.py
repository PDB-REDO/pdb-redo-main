import os
import re
import sys
import json
import numpy as np


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
    
    

def run_dssr(pdb_file, dssr, output_dir):
    path = pdb_file.split("/")
    file_name = path[-1]
    out_name = file_name + ".json"
#    os.chdir(output_dir)
    output_file = output_dir + "/" + out_name
    search_command = (dssr + " -i=" + pdb_file +  " -o=" + output_file + " --json"  + " --more" + " --non-pair > /dev/null 2>&1")
    os.system(search_command)
    return(output_file)

def basepairs(json_file):
    #Nucleiotide types
    DNA = ["DA", "DC", "DG", "DT", "DU"]
    DNAnc = ["08Q","0AD","0DC","0DG","0DT","1AP","1CC","1W5","1WA","2AR","2AT","2BD","2BT","2BU","2DA","2EG","2GT","2NT","2OT","2PR","2ST","3DR","47C","4PC","4PD","4PE","5CM","5FC","5HC","5HU","5IU","5NC","5PC","5PY","5SE","64T","6HA","6HC","6HG","6HT","6MA","6OG","7DA","7GU","85Y","8FG","8MG","8OG","A5L","ABR","ABS","AF2","AFG","ASU","B7C","BGM","BOE","BRU","C34","C6G","C7R","C7S","CAR","CB2","CBR","CDW","CFZ","CTG","D3","D33","DA","DC","DDG","DDN","DFT","DG","DI","DN","DNR","DOC","DRZ","DT","DU","DUZ","DXD","DZM","E","E1X","EDA","EDC","EIT","FDG","FMG","G35","G36","GF2","GMU","GN7","GNE","GSR","GSS","HEU","HOL","IGU","IMC","LCG","LGP","M1G","MA7","MBZ","MG1","MRG","MTR","NCU","NCX","NDU","NMS","NMT","NTT","P","P2T","P2U","PDU","PGN","PPW","PRN","R","S2M","S6G","SOS","T2T","T49","T4S","T5O","T5S","TA3","TAF","TC1","TDY","TFE","TLC","TLN","TTM","U2N","UFT","UMS","US3","USM","X","X4A","XAE","XAL","XCL","XCR","XCS","XCY","XGA","XGL","XGR","XTL","XTR","XTY","XUA","XUG","Y","YRR","Z","ZDU"]
    RNA = ["A", "C", "G", "T", "U"]
    RNAnc = ["0A","0C","0G","0U","0U1","125","126","127","12A","1MA","1MG","23G","2AU","2MA","2MG","2MU","2OM","4AC","4OC","4SU","5BU","5CF","5IC","5MC","5MU","6IA","6NW","6OP","70U","7MG","A","A23","A2L","A2M","A44","A5M","A5O","A6A","A6C","A6G","A6U","AET","AP7","AVC","BMP","C","C2L","C43","CCC","CH","F2T","FHU","G","G25","G2L","G48","G7M","GAO","H2U","I","IC","IG","IU","KAG","LCC","M2G","M5M","MA6","MGV","MIA","MNU","N","N5M","N6G","NF2","O2G","OMC","OMG","OMU","ONE","P5P","PGP","PQ1","PSU","PYO","QUO","RSQ","S4C","SMT","SSU","SUR","T","T6A","U","U2L","U36","U8U","UAR","UOB","UR3","YG","YYG"] 

    #Shear parameters
    mean_shear_at_dna = 0.042
    std_shear_at_dna = 0.095
    mean_shear_at_mix = 0.103
    std_shear_at_mix = 0.060
    
    mean_shear_au_rna = 0.042
    std_shear_au_rna = 0.107
    mean_shear_au_mix = 0.007
    std_shear_au_mix = 0.050  
    
    mean_shear_gc_dna = -0.215
    std_shear_gc_dna = 0.101
    mean_shear_gc_rna = -0.209
    std_shear_gc_rna = 0.115
    mean_shear_gc_mix = -0.175
    std_shear_gc_mix = 0.078
    
    #Stretch parameters
    mean_stretch_at_dna = -0.111
    std_stretch_at_dna = 0.048
    mean_stretch_at_mix = -0.135
    std_stretch_at_mix = 0.027
    
    mean_stretch_au_rna = -0.103
    std_stretch_au_rna = 0.041
    mean_stretch_au_mix = -0.133
    std_stretch_au_mix = 0.020    
    
    mean_stretch_gc_dna = -0.130
    std_stretch_gc_dna = 0.063
    mean_stretch_gc_rna = -0.131
    std_stretch_gc_rna = 0.056
    mean_stretch_gc_mix = -0.150
    std_stretch_gc_mix = 0.037
    
    #Buckle parameters
    mean_buckle_at_dna = 1.325
    std_buckle_at_dna = 7.327
    mean_buckle_at_mix = -5.269
    std_buckle_at_mix = 2.679    
    
    mean_buckle_au_rna = -0.720
    std_buckle_au_rna = 5.951
    mean_buckle_au_mix = 2.380
    std_buckle_au_mix = 2.356
    
    mean_buckle_gc_dna = -0.385
    std_buckle_gc_dna = 8.453
    mean_buckle_gc_rna = -3.403
    std_buckle_gc_rna = 5.437
    mean_buckle_gc_mix = 0.838
    std_buckle_gc_mix = 7.368
    
    #Propeller-twist parameters
    mean_propeller_at_dna = -10.495
    std_propeller_at_dna = 6.315
    mean_propeller_at_mix = -9.742
    std_propeller_at_mix = 3.033    
    
    mean_propeller_au_rna = -11.560
    std_propeller_au_rna = 4.560
    mean_propeller_au_mix = -9.781
    std_propeller_au_mix = 2.039   
    
    mean_propeller_gc_dna = -6.432
    std_propeller_gc_dna = 7.498
    mean_propeller_gc_rna = -10.979
    std_propeller_gc_rna = 5.268
    mean_propeller_gc_mix = -7.512
    std_propeller_gc_mix = 4.648

    lzshear = []
    lzstretch = []
    lzbuckle = []
    lzpropeller = []
    lzbpg = []
    n = 0
    with open(json_file) as file:
        data = json.load(file)
        if "pairs" in data:
            #Print start
            print ("Individual base pair Z-scores")
            print ("-----------------------------")
            print ("Pair                 Shear   Stretch    Buckle Propeller      bpG")
            
            for p in data["pairs"]:
                bp = p["bp"]
                name = p["name"]
                #Only work on WC base pairs
                if name == "WC":
                    try:
                        shear =  float(p["simple_Shear"])
                    except:
                        #Give up on basepair completely
                        continue
                    stretch =  float(p["simple_Stretch"])
                    buckle =  float(p["simple_Buckle"])
                    propeller =  float(p["simple_Propeller"])
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

                    if "a" in bp.lower() and "t" in bp.lower():
                        #Deal with the directionality of shear and buckle
                        if bp.lower()[0:1] == "t":
                            shear = -1*shear
                            buckle = -1*buckle
                        if typ1 == "DNA" and typ2 == "DNA":
                            z_shear = (shear - mean_shear_at_dna) / std_shear_at_dna
                            z_stretch = (stretch - mean_stretch_at_dna) / std_stretch_at_dna
                            z_buckle = (buckle - mean_buckle_at_dna) / std_buckle_at_dna
                            z_propeller = (propeller - mean_propeller_at_dna) / std_propeller_at_dna
                        else:
                            z_shear = (shear - mean_shear_at_mix) / std_shear_at_mix
                            z_stretch = (stretch - mean_stretch_at_mix) / std_stretch_at_mix
                            z_buckle = (buckle - mean_buckle_at_mix) / std_buckle_at_mix
                            z_propeller = (propeller - mean_propeller_at_mix) / std_propeller_at_mix
                    elif "a" in bp.lower() and "u" in bp.lower():
                        #Deal with the directionality of shear and buckle
                        if bp.lower()[0:1] == "u":
                            shear = -1*shear
                            buckle = -1*buckle
                        if typ1 == "RNA" and typ2 == "RNA":
                            z_shear = (shear - mean_shear_au_rna) / std_shear_au_rna
                            z_stretch = (stretch - mean_stretch_au_rna) / std_stretch_au_rna
                            z_buckle = (buckle - mean_buckle_au_rna) / std_buckle_au_rna
                            z_propeller = (propeller - mean_propeller_au_rna) / std_propeller_au_rna  
                        else:
                            z_shear = (shear - mean_shear_au_mix) / std_shear_au_mix
                            z_stretch = (stretch - mean_stretch_au_mix) / std_stretch_au_mix
                            z_buckle = (buckle - mean_buckle_au_mix) / std_buckle_au_mix
                            z_propeller = (propeller - mean_propeller_au_mix) / std_propeller_au_mix
                    elif "g" in bp.lower() and "c" in bp.lower():
                        #Deal with the directionality of shear and buckle
                        if bp.lower()[0:1] == "c":
                            shear = -1*shear
                            buckle = -1*buckle
                        if typ1 == "DNA" and typ2 == "DNA":
                            z_shear = (shear - mean_shear_gc_dna) / std_shear_gc_dna
                            z_stretch = (stretch - mean_stretch_gc_dna) / std_stretch_gc_dna
                            z_buckle = (buckle - mean_buckle_gc_dna) / std_buckle_gc_dna
                            z_propeller = (propeller - mean_propeller_gc_dna) / std_propeller_gc_dna   
                        elif typ1 == "RNA" and typ2 == "RNA":
                            z_shear = (shear - mean_shear_gc_rna) / std_shear_gc_rna
                            z_stretch = (stretch - mean_stretch_gc_rna) / std_stretch_gc_rna
                            z_buckle = (buckle - mean_buckle_gc_rna) / std_buckle_gc_rna
                            z_propeller = (propeller - mean_propeller_gc_rna) / std_propeller_gc_rna  
                        else:
                            z_shear = (shear - mean_shear_gc_mix) / std_shear_gc_mix
                            z_stretch = (stretch - mean_stretch_gc_mix) / std_stretch_gc_mix
                            z_buckle = (buckle - mean_buckle_gc_mix) / std_buckle_gc_mix
                            z_propeller = (propeller - mean_propeller_gc_mix) / std_propeller_gc_mix                          

                    #Calculate ZbgG
                    z_bpg = (abs(z_shear)+abs(z_stretch)+abs(z_buckle)+abs(z_propeller))/4
                    
                    #Fill the lists
                    lzshear.append(z_shear)
                    lzstretch.append(z_stretch)
                    lzbuckle.append(z_buckle)
                    lzpropeller.append(z_propeller)
                    lzbpg.append(z_bpg)
                    
                    #Print individual scores
                    print ("{:16s}".format(p["nt1"]+ "-" + p["nt2"])+ "{:>10.2f}{:>10.2f}{:>10.2f}{:>10.2f}{:>10.2f}".format(z_shear,z_stretch,z_buckle,z_propeller,z_bpg)) 

                    n += 1
        else:
            print (" ")
            print ("Overall base pair rmsZ-scores")
            print ("-----------------------------")  
            print ("Shear rmsZ    : NA") 
            print ("Stretch rmsZ  : NA")
            print ("Buckle rmsZ   : NA")
            print ("Propeller rmsZ: NA")
            print ("bpG rmsZ      : NA")
            print ("Count         : NA")
            return ""

        #Cast the data
        azshear = np.array(lzshear)
        azstretch = np.array(lzstretch)
        azbuckle = np.array(lzbuckle)
        azpropeller = np.array(lzpropeller)
        azbpg = np.array(lzbpg)

        print (" ")
        print ("Overall base pair rmsZ-scores (sigma)")
        print ("-----------------------------")       
        print ("Shear rmsZ    : " +  "{:.3f}".format(rms(azshear)) + " (" + "{:.3f}".format(jackknife_stdev(azshear, rms)) +")") 
        print ("Stretch rmsZ  : " +  "{:.3f}".format(rms(azstretch)) + " (" + "{:.3f}".format( jackknife_stdev(azstretch, rms)) +")")
        print ("Buckle rmsZ   : " +  "{:.3f}".format(rms(azbuckle)) + " (" + "{:.3f}".format( jackknife_stdev(azbuckle, rms)) +")")
        print ("Propeller rmsZ: " +  "{:.3f}".format(rms(azpropeller)) + " (" + "{:.3f}".format( jackknife_stdev(azpropeller, rms)) +")")
        print ("bpG rmsZ      : " +  "{:.3f}".format(rms(azbpg)) + " (" + "{:.3f}".format( jackknife_stdev(azbpg, rms)) +")")
        print ("Count: ",n)
        return ""

#Get input file
pdb_file = os.path.abspath(sys.argv[1])
dssr_path = os.path.abspath(sys.argv[2])
search_command = ('mkdir -p dssr_output')
os.system(search_command)
dssr_file = run_dssr(pdb_file, dssr_path, 'dssr_output')
rmsz_values = basepairs(dssr_file)
print(rmsz_values)
