from logging import ERROR, error
from pathlib import Path

import json
import gzip 
import re 
import sys

class ligandvalidation:
    def __init__(self, input_path):
        self.input_path = input_path

##### Set  up file system and files required #####
    
    def define_pdbid(versions_file):
        try:
            versions_file_open = Path(versions_file).read_text() 
            versions = json.loads(versions_file_open)  
            pdb_id = versions['data']['PDBID']
        except:
            sys.exit('Error while reading PDBID from versions.json')
        return pdb_id


    def define_path_files(self, input_path): 
        '''Set all paths and files needed, and check if all exist'''

        work_path = Path(input_path)
        if work_path.is_dir() == False:
            sys.exit('Work_path (tempdir) could not be found, stopped running script') 
        else:
            work_path = work_path

        ligval_files = work_path.glob('ligval_*.log')
        ligval_list = list(ligval_files)
        if ligval_list == []:
            sys.exit('No ligand validation files present, stopped running script as no ligands reported for this structure model')
        else:
            ligval_list = ligval_list 

        versions_file = work_path/'versions.json'
        if versions_file.exists() == False:
            sys.exit('versions_file could not be found, stopped running script')
        else:    
            pdb_id = self.define_pdbid(versions_file)


        jsonfile_0cyc = work_path/f'{pdb_id}_0cyc.json'
        if jsonfile_0cyc.exists() == False:
            sys.exit('jsonfile_0cyc could not be found, stopped running script')
        else: 
            jsonfile_0cyc = jsonfile_0cyc 


        jsonfile_final = work_path/f'{pdb_id}_final.json'
        if jsonfile_final.exists() == False:
            sys.exit('jsonfile_final could not be found, stopped running script')
        else: 
            jsonfile_final = jsonfile_final 

        pdbfile_final = work_path/f'{pdb_id}_final.pdb'
        if pdbfile_final.exists() == False:
            sys.exit('pdbfile_final could not be found, stopped running script')
        else: 
            pdbfile_final = pdbfile_final

        output_file = work_path/f'{pdb_id}_ligval.json'


        return ligval_list, versions_file, jsonfile_0cyc, jsonfile_final, output_file, pdbfile_final

##### All functions required to obtain the "Data"-block of the output file #####
    def get_data(self, versions_file):
        '''Get all the data for the "Data"-block at top of the PDBID_ligval.json file, return as dictionary'''
        ### Open and read versions.json ###
        versions_file_open = Path(versions_file).read_text() 
        versions = json.loads(versions_file_open)    
        
        ### Create dictionary containing all extracted above and return this dictionary ###
        try:
            dict_data = { "pdb_id" : versions['data']['PDBID'], 
                          "coordinates_revision_date_pdb" : versions['data']['coordinates_revision_date_pdb'], 
                          "coordinates_revision_major_mmCIF" : versions['data']['coordinates_revision_major_mmCIF'], 
                          "coordinates_revision_minor_mmCIF" : versions['data']['coordinates_revision_minor_mmCIF'], 
                          "coordinates_edited" : versions['data']['coordinates_edited'],
                          "reflections_revision" : versions['data']['reflections_revision'], 
                          "reflections_edited" : versions['data']['reflections_edited'], 
                          "pdb_redo_version" : versions['software']['pdb-redo']['version'], 
                          "pdb_redo_used" : versions['software']['pdb-redo']['used'] } 
        except KeyError:
            raise sys.exit('KeyError while reading from versions.json, stopped running script')
        return dict_data

##### All functions required to obtain the "Ligand validation data"-block of the output file #####
    def get_pdb_details_ligand(self, jsonfile_final, ligand):
        ligand_open = Path(ligand).read_text()
        residue = ligand_open.split('Residue')
        pdb_compID = residue[1].split()[0]
        pdb_seqNum = int(residue[1].split()[2])
        pdb_strandID = residue[1].split()[1]

        jsonfile_final_open = Path(jsonfile_final).read_text() 
        final = json.loads(jsonfile_final_open)
        for element in final:
            if element['pdb']['compID'] == pdb_compID and element['pdb']['seqNum'] == pdb_seqNum and element['pdb']['strandID'] == pdb_strandID:
                asymID_mmcif = element['asymID']
                compID_mmcif = element['compID']
                seqID_mmcif = element['seqID']
                pdb_insCode = element['pdb']['insCode']
        
        return asymID_mmcif, compID_mmcif, seqID_mmcif, pdb_compID, pdb_insCode, pdb_seqNum, pdb_strandID    

    def was_ligand_added(self, ligand):
        ligand_open = Path(ligand).read_text()
        if 'ori' in ligand_open:
            return False            ### Ligand does exist in original PDB entry ###
        else:
            return True             ### Ligand was added by PDB-REDO ###

    def get_added_ligand_occupancy(self, pdb_compID, pdb_seqNum, pdb_strandID, pdbfile_final):
        ''' Only when ligand was added by PDB-REDO '''
        pdbfile_open = Path(pdbfile_final).read_text().split('\n')
        atoms_list = [line for line in pdbfile_open if line[0:5] == 'ATOM' or line[0:6] == 'HETATM']
        occupancy_list = [atom[54:60] for atom in atoms_list if len(atom) == 80 and pdb_compID == atom[17:20] and str(pdb_seqNum) in atom[22:26] and pdb_strandID == atom[21]]
        if occupancy_list != []:
            occupancy = float(occupancy_list[0])
        else:
            sys.exit(f'Something went wrong in getting occupancy for {pdb_compID, pdb_seqNum, pdb_strandID}. Most likely something wrong with the input or the ligand not present in pdbfile_final')
        return occupancy
        


    def get_rmsd_shift(self, ligand):
        ''' Only when ligand was not added by PDB-REDO '''     
        try:   
            all_atom_rmsd = [float(line.split()[10]) for line in Path(ligand).read_text().split('\n') if re.match('Residue', line) and re.search(' A RMSD over ', line)][0]
            
            try:
                shifted_list = [float(line.split()[14]) for line in Path(ligand).read_text().split('\n') if re.match('Atom ', line) and re.search(' A RMSD over ', line)]
                shifted = sum([1 for shift in shifted_list if shift > 0.5])
            except IndexError:
                shifted = None      ### Error in either script or YASARA logfile ###
        
        except IndexError:          ### No rmsd reported in YASARA logfile, probabaly due to error in YASARA ###
            all_atom_rmsd = None
            shifted = None
        return all_atom_rmsd, shifted
    
    
    def get_density_fit_ori(self, jsonfile_final, jsonfile_0cyc, ligand):
        asymID_mmcif, compID_mmcif, seqID_mmcif, pdb_compID, pdb_insCode, pdb_seqNum, pdb_strandID = self.get_pdb_details_ligand(self, jsonfile_final, ligand)
        try:
            jsonfile_0cyc_open = Path(jsonfile_0cyc).read_text() 
            file_0cyc = json.loads(jsonfile_0cyc_open)
        except:
            sys.exit('error while reading and loading jsonfile_0cyc, stopped running script')

        for element in file_0cyc:
            if element['pdb']['compID'] == pdb_compID and element['pdb']['seqNum'] == pdb_seqNum and element['pdb']['strandID'] == pdb_strandID:
                RSRfactor = float(round(element['RSR'], 3))
                RSCCS = float(round(element['RSCCS'], 3))
                EDIAm = float(round(element['EDIAm'], 3))
                OPIA = float(round(element['OPIA'], 3))
                
        density_fit = {}
        density_fit['real_space_Rfactor'] = RSRfactor
        density_fit['real_space_correlation'] = RSCCS
        density_fit['EDIAm_density_fit'] = EDIAm
        density_fit['OPIA_density_coverage'] = OPIA
        return density_fit

    def get_density_fit_redo(self, jsonfile_final, ligand):
        asymID_mmcif, compID_mmcif, seqID_mmcif, pdb_compID, pdb_insCode, pdb_seqNum, pdb_strandID = self.get_pdb_details_ligand(self, jsonfile_final, ligand)
        try:
            jsonfile_final_open = Path(jsonfile_final).read_text() 
            final = json.loads(jsonfile_final_open)
        except:
            sys.exit('error while reading and loading jsonfile_final, stopped running script')

        for element in final:
            if element['pdb']['compID'] == pdb_compID and element['pdb']['seqNum'] == pdb_seqNum and element['pdb']['strandID'] == pdb_strandID:
                RSRfactor = float(round(element['RSR'], 3))
                RSCCS = float(round(element['RSCCS'], 3))
                EDIAm = float(round(element['EDIAm'], 3))
                OPIA = float(round(element['OPIA'], 3))
                
        density_fit = {}
        density_fit['real_space_Rfactor'] = RSRfactor
        density_fit['real_space_correlation'] = RSCCS
        density_fit['EDIAm_density_fit'] = EDIAm
        density_fit['OPIA_density_coverage'] = OPIA
        return density_fit

    def get_heat_of_formation(self, ligand, var):                  
        ''' var = 1 ; for original model or ligand added by PDB-REDO '''
        ''' var = 2 ; for PDB-REDO model when ligand was not added by PDB-REDO '''
        try:
            energy = [float(line.split()[7]) for line in Path(ligand).read_text().split('\n') if re.match(f'Energy of formation in object {var}', line)][0] 

        except IndexError:         ### No HOF reported in YASARA logfile, probabaly due to error in YASARA ###
            energy = None
        
        heat_of_formation = {}
        heat_of_formation['energy'] = energy
        heat_of_formation['energy_unit'] = 'kJ/mol'
        return heat_of_formation

    def get_interactions(self, ligand, var):
        ''' var = ori ; for original model  '''
        ''' var = new ; for PDB-REDO model, also when ligand was added by PDB-REDO '''
        ### Functions to get all the interactions ###        
        def get_bumps(interaction):
            try:
                bumps = sum([float(item.split()[-1]) for item in interaction if 'Contacts to Residue' in item])
            except ValueError:
                bumps = 0
            return bumps
                        
        def get_hbond(interaction):
            try:
                hbond_count = float(interaction[-3].split()[0])
                try:
                    hbond_strength = float(interaction[-3].split()[-2]) * -1 

                except ValueError:
                    hbond_strength = None
                    raise f'Error in extracting hbond_strength for {ligand}'
            except ValueError:
                hbond_count = 0
                hbond_strength = None
            except IndexError:
                hbond_count = 0
                hbond_strength = None
            return hbond_count, hbond_strength

        def get_hydrophobic(interaction):
            try: 
                hydpho_count = sum([float(item.split(':')[1].split()[0]) for item in interaction if 'interactions with strength' in item])
                try:
                    hydpho_strength = float(round( sum([float(item.split(':')[1].split()[-1]) for item in interaction if 'interactions with strength' in item]) , 3))
                except ValueError:
                    hydpho_strength = None
                    raise f'Error in calculation hydpho_strength for {ligand}'
            except ValueError:
                hydpho_count = 0
                hydpho_strength = None
            return hydpho_count, hydpho_strength

        def get_pipi(interaction):
            try:
                pipi_count = sum([float(item.split(':')[1].split()[0]) for item in interaction if 'interactions with strength' in item])
                if pipi_count != 0 :
                    try:
                        pipi_strength = float(round( sum([float(item.split(':')[1].split()[-1]) for item in interaction if 'interactions with strength' in item]), 3))  
                    except ValueError:
                        pipi_strength = None
                        raise f'Error in calculation pipi_strength for {ligand}'
                else:
                    pipi_strength = None

            except ValueError:
                pipi_count = 0
                pipi_strength = None
            return pipi_count, pipi_strength

        def get_catpi(interaction):
            try:
                catpi_count = sum([float(item.split(':')[1].split()[0]) for item in interaction if 'interactions with strength' in item])
                if catpi_count != 0 :
                    try:
                        catpi_strength = float(round( sum([float(item.split(':')[1].split()[-1]) for item in interaction if 'interactions with strength' in item]) , 3))
                    except ValueError:
                        catpi_strength = None
                        raise f'Error in calculation catpi_strength for {ligand}'
                else:
                    catpi_strength = None
            except ValueError:
                catpi_count = 0
                catpi_strength = None
            return catpi_count, catpi_strength

        ### Get all the counts and strengths of the interactions that involve the ligand ###
        ligand_log = Path(ligand).read_text().split('Start')
        bumps = get_bumps([interaction.split('\n') for interaction in ligand_log if f'bumps {var}' and f'End bumps {var}' in interaction][0])
        hbond_count, hbond_strength = get_hbond([interaction.split('\n') for interaction in ligand_log if f'hbond {var}' and f'End hbond {var}' in interaction][0])
        hydpho_count, hydpho_strength = get_hydrophobic([interaction.split('\n') for interaction in ligand_log if f'hydpho {var}' and f'End hydpho {var}' in interaction][0])
        pipi_count, pipi_strength = get_pipi([interaction.split('\n') for interaction in ligand_log if f'pipi {var}' and f'End pipi {var}' in interaction][0])
        catpi_count, catpi_strength = get_catpi([interaction.split('\n') for interaction in ligand_log if f'catpi {var}' and f'End catpi {var}' in interaction][0])

        ### Make dictionary containing interaction details for ligand ori ###
        interactions = {}
        interactions['bumps_count'] = bumps
        interactions['H_bonds_count'] = hbond_count
        interactions['H_bonds_strength'] = hbond_strength
        interactions['H_bonds_unit'] = 'kJ/mol'
        interactions['hydrophobic_count'] = hydpho_count
        interactions['hydrophobic_strength'] = hydpho_strength
        interactions['hydrophobic_unit'] = 'knowledge_based_potential'
        interactions['pi_pi_count'] = pipi_count
        interactions['pi_pi_strength'] = pipi_strength
        interactions['pi_pi_unit'] = 'knowledge_based_potential'
        interactions['cation_pi_count'] = catpi_count
        interactions['cation_pi_strength'] = catpi_strength
        interactions['cation_pi_unit'] = 'knowledge_based_potential'
        
        return interactions
      
 
    def get_ligand_validation(self, ligval_list, jsonfile_final, jsonfile_0cyc, pdbfile_final):
        '''Get all the data for all ligands present for the "Ligand_validation_data"-block in the PDBID_ligval.json file'''      
        all_ligands_validation = []
        for ligand in ligval_list:
            ### Get ligand properties and add to dictionary ###
            asymID_mmcif, compID_mmcif, seqID_mmcif, pdb_compID, pdb_insCode, pdb_seqNum, pdb_strandID = self.get_pdb_details_ligand(self, jsonfile_final, ligand)

            ligand_validation = {}
            ligand_validation['asymID'] = asymID_mmcif
            ligand_validation['compID'] = compID_mmcif
            ligand_validation['seqID'] = seqID_mmcif
            ligand_validation['pdb'] = {}
            ligand_validation['pdb']['compID'] = pdb_compID
            ligand_validation['pdb']['insCode'] = pdb_insCode
            ligand_validation['pdb']['seqNum'] = pdb_seqNum
            ligand_validation['pdb']['strandID'] = pdb_strandID
            
            ### Define if ligand was added by PDB-REDO or not ###
            ligand_validation['added'] = self.was_ligand_added(self, ligand)
            
            if ligand_validation['added'] == False:
            ### Get data from original model and PDB-REDO model if added = False ### 
                all_atom_rmsd, shifted = self.get_rmsd_shift(self, ligand)

                ligand_validation['all_atom_rmsd_in_A'] = all_atom_rmsd
                ligand_validation['atoms_shifted_more_than_0.5A'] = shifted
                ligand_validation['original_model'] = {}
                ligand_validation['original_model']['density_fit'] = self.get_density_fit_ori(self, jsonfile_final, jsonfile_0cyc, ligand)
                ligand_validation['original_model']['heat_of_formation'] = self.get_heat_of_formation(self, ligand, 1)
                ligand_validation['original_model']['interactions'] = self.get_interactions(self, ligand, 'ori') 
                ligand_validation['pdb_redo_model'] = {}
                ligand_validation['pdb_redo_model']['density_fit'] = self.get_density_fit_redo(self, jsonfile_final, ligand)
                ligand_validation['pdb_redo_model']['heat_of_formation'] = self.get_heat_of_formation(self, ligand, 2)
                ligand_validation['pdb_redo_model']['interactions'] = self.get_interactions(self, ligand, 'new')     
            else:              
            ### Get data from PDB-REDO model (if added = True ) ###
                ligand_validation['pdb_redo_model'] = {}
                ligand_validation['pdb_redo_model']['refined_occupancy'] = self.get_added_ligand_occupancy(self, pdb_compID, pdb_seqNum, pdb_strandID, pdbfile_final)
                ligand_validation['pdb_redo_model']['density_fit'] = self.get_density_fit_redo(self, jsonfile_final, ligand)
                ligand_validation['pdb_redo_model']['heat_of_formation'] = self.get_heat_of_formation(self, ligand, 1)
                ligand_validation['pdb_redo_model']['interactions'] = self.get_interactions(self, ligand, 'new')
            
            all_ligands_validation.append(ligand_validation)
        return all_ligands_validation


##### Combine data and ligand_validation_data blocks in one json formatted outputfile #####
    def main(self, input_path):
        '''Generate output file PDBID_ligval.json containing all ligand validation data'''
        ligval_list, versions_file, jsonfile_0cyc, jsonfile_final, output_file, pdbfile_final = self.define_path_files(self, input_path)
        pdb_id = self.define_pdbid(versions_file)
        data = self.get_data(self, versions_file)
        ligand_validation = self.get_ligand_validation(self, ligval_list, jsonfile_final, jsonfile_0cyc, pdbfile_final)
        
        combined = {}
        combined['Data'] = data
        combined['Ligand_validation_data'] = ligand_validation
        with open(output_file, 'w') as f:
            json.dump(combined, f, indent = 2, sort_keys = False)        
        
        print(f'Generated {pdb_id}_ligval.json')
        return output_file


##### Initiate running script #####
if __name__ == '__main__':
        input_path = sys.argv[1]
        ligandvalidation.main(ligandvalidation, input_path)



