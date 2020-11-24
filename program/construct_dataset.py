#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:22:20 2020

@author: kacper
"""
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB import PDBList
from Bio.PDB.DSSP import DSSP


#%%
def ID_to_URL (ID): #ticked
    """Input  : Prot ID
   Return : Url to *.pdb file
   """
    pdb_id = "https://files.rcsb.org/download/%s.pdb"
    return (pdb_id % ID.rstrip())

#%%

def generate_model_and_dssp_object(protein_name):#ticked
    """
    

    Parameters
    ----------
    protein_name : str
        protein ID, to acces PDB database

    Returns
    -------
    structure : Bio.PDB.Structure.Structure
        The object representing a model in a structure.
        The Structure class contains a collection of Model instances.
    dssp_object : Bio.PDB.DSSP.DSSP
        parsed secondary structure and accessibility.
    """

    pdbl = PDBList()
    #Quick access to the structure lists on the PDB or its mirrors. 
    
    pdb_filepath = pdbl.retrieve_pdb_file(protein_name, file_format='pdb')
    #Fetch PDB structure file from PDB server, and store it locally.
    
    parser = PDBParser()
    #Parser for PDB files. read docs, for class members. usefu;
    
    structure = parser.get_structure(protein_name, pdb_filepath)
    # returns a structure 
    
    model = structure[0]
    #to get Model id 
    
    dssp_object = make_a_dssp_model(model, pdb_filepath)
    #constructs an dssp object
    
    return structure, dssp_object

def make_a_dssp_model(model, pdb_filepath, dssp='dssp'):#ticked
    """
    
    Parameters
    ----------
    model : TYPE
        The Structure class contains a collection of Model instances.
    pdb_filepath : str
        DESCRIPTION.
    dssp : str, optional
        DESCRIPTION. The default is 'dssp'.

    Returns
    -------
    dssp_object : Bio.PDB.DSSP.DSSP
        parsed secondary structure and accessibility.

    """
    #DSSP class, which maps Residue objects to their secondary structure
    try:
        dssp_object = DSSP(model, pdb_filepath, dssp)
    except Exception as e:
        if type(e) is not Exception:
            raise
        print('oops')
        
    
    return dssp_object
#%%
def gen_array_of_CA_coords(structure, chain_id=None): #ticked
    """
    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The Structure class contains a collection of Model instances.
        The object representing a model in a structure.

    Returns
    -------
    array_of_coords, numpy.array
        returns a flatten numpy.array, composed of coordinates of
        alpha carbons, formatted using PCA principle, as follows
        [x_1 y_1 z_1 x_2 y_2 z_2 ... x_n y_n z_n]

    """
    if chain_id:
        array_of_coords = np.array([atom.get_coord() for atom \
                                in structure[0][chain_id].get_atoms()\
                                if atom.get_name() == "CA" ])
    else:    
        array_of_coords = np.array([atom.get_coord() for atom \
                                in structure[0].get_atoms()\
                                if atom.get_name() == "CA" ])
        
    return array_of_coords.flatten()
#%%
def make_onehot_encoding_of_secondary_structure_as_q8(dssp_model, chain_id=None): #ticked
    """
    Parameters
    ----------
    dssp_model : Bio.PDB.DSSP.DSSP object

    Returns
    -------
    as Q8
    onehot_encoded_secondary_structure : TYPE == numpy.array
        Numpy array of one hot encoded secondary structures

    """
    
    """
    H == Alpha helix (4-12)    B == Isolated beta-bridge residue
    E == Strand                G == 3-10 helix
    I == Pi helix              T == Turn
    S == Bend                  - == None
    """
    mapping = {'H':0, 'B':1, 'E':2, 'G':3, 'I':4, 'T':5, 'S':6, '-':7}
    """
    mapping[indeq] for indeq in x will change ss index to given number
    [dssp[key][2] for key in list(dssp.keys())] will return a list
    of ss indexes, as ss value is at index 2 of dssp tuple
    last make that list of mapped values to numpy array
    and then one-hot encode it
    """
    if chain_id:
        list_of_keys_in_chain = [key for key in list(dssp_model.keys()) \
                             if key[:][0] == chain_id ]
            
        onehot_encoded_ss = onehot(np.array( [ mapping[indeq] for indeq in \
                    [dssp_model[key][2] for key in list_of_keys_in_chain]]))
    else:        
        #print([ mapping[indeq] for indeq in \
         #           [dssp_model[key][2] for key in list(dssp_model.keys())]])
        onehot_encoded_ss = onehot(np.array( [ mapping[indeq] for indeq in \
                    [dssp_model[key][2] for key in list(dssp_model.keys())]]))
    """
    generates list of keys to acces dssp_object_attributes, as follows :
    (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
    NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
    NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
    """
    return onehot_encoded_ss

def make_onehot_encoding_of_secondary_structure_as_Q3(dssp_model, chain_id=None):
    mapping = {'H':0, 'G':0, 'I':0, 'E':1, 'B':1, 'T':2, 'S':2, '-':2}
    
    dssp_keys = list(dssp_model.keys())
    if chain_id:
        list_of_keys_in_chain = [key for key in dssp_keys \
                                 if key[:][0] == chain_id ]
        
        onehot_encoded_ss_q3 = onehot(np.array( [ mapping[index] for index in \
                [dssp_model[key][2] for key in list_of_keys_in_chain]]))
                                                 
    else:
        onehot_encoded_ss_q3 = onehot(np.array( [ mapping[indeq] for indeq in \
                    [dssp_model[key][2] for key in dssp_keys]]))

    
    return onehot_encoded_ss_q3

def onehot(arr): #ticked
    """
    Parameters
    ----------
    arr : numpy.array
        numpy array of secondary structure affiliation, encoded in numbers
        varying from 0 to 7, stated in <make_onehot_encoding> function

    Returns
    -------
    out : numpy.array
        onehot encoded numpy array, of 2 dimensions, 8 columns each for ss
        , and as many rows as there are records in input arr

    """
    out = np.zeros((arr.shape[0], 8))
    for i, k in enumerate(arr):
        out[i, k-1] = 1
    return out

def make_regular_encoding_of_secondary_structure_as_Q8(dssp_model, chain_id=None):#ticked
    """
    H == Alpha helix (4-12)    B == Isolated beta-bridge residue
    E == Strand                G == 3-10 helix
    I == Pi helix              T == Turn
    S == Bend                  - == None
    """
    mappingQ8 = {'H':0, 'B':1, 'E':2, 'G':3, 'I':4, 'T':5, 'S':6, '-':7}
    """
    mapping[indeq] for indeq in x will change ss index to given number
    [dssp[key][2] for key in list(dssp.keys())] will return a list
    of ss indexes, as ss value is at index 2 of dssp tuple
    last make that list of mapped values to numpy array
    and then one-hot encode it
    """

    """
    dssp_model.keys() - > generates list of keys 
    to acces dssp_object_attributes, as follows :
    (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
    NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
    NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
    """
    if chain_id:
        list_of_keys_in_chain = [key for key in list(dssp_model.keys()) \
                             if key[:][0] == chain_id ]

        regular_encoded_ssq8 = np.array( [ mappingQ8[indeq] for indeq in \
                    [dssp_model[key][2] for key in list_of_keys_in_chain]])
     
    else:       
        regular_encoded_ssq8 = np.array( [ mappingQ8[indeq] for indeq in \
                    [dssp_model[key][2] for key in list(dssp_model.keys())]])
    
    
    return regular_encoded_ssq8


def make_regular_encoding_of_secondary_structure_as_Q3(dssp_model, chain_id=None):#ticked
    """
    H == Alpha helix (4-12)    B == Isolated beta-bridge residue
    E == Strand                G == 3-10 helix
    I == Pi helix              T == Turn
    S == Bend                  - == None
    """
    mappingQ3 = {'H':0, 'G':0, 'I':0, 'E':1, 'B':1, 'T':2, 'S':2, '-':2}
    """
    mapping[indeq] for indeq in x will change ss index to given number
    [dssp[key][2] for key in list(dssp.keys())] will return a list
    of ss indexes, as ss value is at index 2 of dssp tuple
    last make that list of mapped values to numpy array
    and then one-hot encode it
    """

    """
    dssp_model.keys() - > generates list of keys 
    to acces dssp_object_attributes, as follows :
    (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
    NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
    NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
    """
    if chain_id:
        list_of_keys_in_chain = [key for key in list(dssp_model.keys()) \
                             if key[:][0] == chain_id ]

        regular_encoded_ssq3 = np.array( [ mappingQ3[indeq] for indeq in \
                    [dssp_model[key][2] for key in list_of_keys_in_chain]])
     
    else:       
        regular_encoded_ssq3 = np.array( [ mappingQ3[indeq] for indeq in \
                    [dssp_model[key][2] for key in list(dssp_model.keys())]])
    
    
    return regular_encoded_ssq3

#%%


def GAME(protein_id, q_type=0,ss_code_type=0):#ticked
    """
    

    Parameters
    ----------
    protein_id : string
        protein id, formatted as id, or or followed by _X,
        where X is chain_id
    q_type : int
        0 for Q8 format, 
        1 for Q3 format
    ss_code_type : int
        0 for regular format
        1 for onehot encoded format

    Returns
    -------
    array_of_CA_as_PCA : numpy array
        array of CA coordinated in PCA format
    ss_as_onehot : numpy array
        array of secondary structure codes , as one-hot encoding

    """
    
    protein_pdb_id = None
    chain_id       = None
    if len(protein_id) == 6:
        chain_id       = protein_id[5]
        protein_pdb_id = protein_id[:4]
    else:
        protein_pdb_id = protein_id
    
    model, dssp = generate_model_and_dssp_object(protein_pdb_id)
    
    
    array_of_CA_as_PCA = None
    
    if chain_id:
        array_of_CA_as_PCA = gen_array_of_CA_coords(model, chain_id)
        
        if q_type == 0 :#for q type
            if ss_code_type == 0:#for ss code type 0 is regular
                ss_encoded = make_regular_encoding_of_secondary_structure_as_Q8(dssp, chain_id)
            else: #else onehot encoded
                ss_encoded = make_onehot_encoding_of_secondary_structure_as_q8(dssp, chain_id)
        
        else:
            if ss_code_type == 0:
                ss_encoded = make_regular_encoding_of_secondary_structure_as_Q3(dssp, chain_id)
            else:
                ss_encoded = make_onehot_encoding_of_secondary_structure_as_Q3(dssp, chain_id)
    
    else:
        array_of_CA_as_PCA = gen_array_of_CA_coords(model)
        
        if q_type == 0:
            if ss_code_type == 0:
                ss_encoded = make_regular_encoding_of_secondary_structure_as_Q8(dssp)
            else:    
                ss_encoded = make_onehot_encoding_of_secondary_structure_as_q8(dssp)
        
        else:           
            if ss_code_type == 0:
                ss_encoded = make_regular_encoding_of_secondary_structure_as_Q8(dssp)
            else:
                ss_encoded = make_onehot_encoding_of_secondary_structure_as_Q3(dssp)
    
    return array_of_CA_as_PCA, ss_encoded



#%%
def gen_dataset(target_list_filepath, q_type=0, ss_code_type=0):#ticked
    """
    

    Parameters
    ----------
    target_list_filepath : str
        path to a file containing entries
    q_type : int
        0 for Q8 format, 
        1 for Q3 format
    ss_code_type : int
        0 for regular format
        1 for onehot encoded format

    Returns
    -------
    A tab for all entries in target_list
    with ca coords as temp_tab[0] and ss_codes as temp_tab[1].

    """    
    
    with open (target_list_filepath) as target_in:
        tab_of_targets = [name.strip('\n').strip() for name in target_in]
    print(tab_of_targets)
        
    temp_tab = [ [], [] ]
    
    
    #dodac ss_code_Type
    for target in tab_of_targets:
        try:
            if q_type==0:
                if ss_code_type == 0:
                    coords, ss_enc = GAME(target)
                else:
                    coords, ss_enc = GAME(target, ss_code_type=1)
            else:
                if ss_code_type == 0:
                    coords, ss_enc = GAME(target, q_type=1)
                else:
                    coords, ss_enc = GAME(target, q_type=1,ss_code_type=1)
            
            temp_tab[0].append(coords)
            temp_tab[1].append(ss_enc)
        except FileNotFoundError:
            pass
    
    return np.array(temp_tab)
    



#test nowego dendata
#%%
fp = "/home/kacper/Dokumenty/calfa/kod/test_names.txt"
datatry_reg_q8 = gen_dataset(fp)
datatry_onh_q8 = gen_dataset(fp, ss_code_type=1)
datatry_reg_q3 = gen_dataset(fp, q_type=1)
datatry_onh_q3 = gen_dataset(fp, q_type=1, ss_code_type=1)


#%%

temp_path = "/home/kacper/Dokumenty/calfa/kod/test_names.txt"

target_path = "/home/kacper/Dokumenty/calfa/kod/protein_list.txt"
main_path = "/home/kacper/Dokumenty/calfa"
datasets_path = "/home/kacper/Dokumenty/calfa/datasets"


Q8_regular_ss_codes_path  = f"{datasets_path}/q8_pcacc_regular_ss.npy"
Q8_onehot_ss_codes_path   = f"{datasets_path}/q8_pcacc_onehot_ss.npy"
Q3_regular_ss_codes_path  = f"{datasets_path}/q3_pcacc_regular_ss.npy"
Q3_onehot_ss_codes_path   = f"{datasets_path}/q3_pcacc_onehot_ss.npy"


np.save(Q8_regular_ss_codes_path,datatry_reg_q8)

np.save(Q8_onehot_ss_codes_path,datatry_onh_q8)

np.save(Q3_regular_ss_codes_path,datatry_reg_q3)

np.save(Q3_onehot_ss_codes_path,datatry_onh_q3)

#%%








    

