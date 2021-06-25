import pandas as pd
import numpy as np
from rdkit import Chem
import multiprocessing as mp
from functools import partial
import logging

log = logging.getLogger("output_log")
log.info("Starting up the log!")

# Load in a sample dataset of Chem SMILES and their HBA/HBD counts
chem_df = pd.read_csv('molecules.csv')
molecules = chem_df['Canonical_SMILES'].unique()

# get a mapping of molecule by index to track its progress in the parallel
# processing shell output
molecule_ix_dict = {}
for molecule_ix in range(len(molecules)):
    molecule_ix_dict[molecules[molecule_ix]] = molecule_ix
# print(chem_df.head())

# Create dataframe that will contain polymerization information
df = pd.DataFrame(columns = ['Chem','HBA','HBD','amide','amide_HBA','amide_HBD','ester','ester_HBA',
                             'ester_HDB'])


# Merge Canonical_SMILES, HBA, and HBD data from hydrogen bond dataset to polymerization dataset
df['Chem'] = pd.Series(chem_df['Canonical_SMILES'])
df['HBA'] = pd.Series(chem_df['HBA'])
df['HBD'] = pd.Series(chem_df['HBD'])
df.head()



# Check how a SMILES string is written in Canonical form. 
mol = Chem.MolToSmiles(Chem.MolFromSmiles('O=COC'))
print(mol)
Chem.MolFromSmiles(mol)

# The function 'submatches' lets us know for a given molecule (mol) how many of a certain substructure (sub) there are.
def submatches(mol,sub):
    substruct_matches = Chem.MolFromSmiles(mol).GetSubstructMatches(Chem.MolFromSmiles(sub))
    # debug CO structure column... CO should not be double-counting stuff
    # if sub == "CO":
    #     print(f"Molecule: {mol}, substructure pattern: {sub}")
    #     print("Substructure matches: ", substruct_matches)
    #     wait = input("Press enter to continue...")
    return substruct_matches


# Specifying substructures to search for in each small molecule. Key below. Each substructure is written in 
# canonical SMILES form, which helps us avoid graph isomorphism (i.e. functional groups are only written in one way in 
# canonical SMILES form)

################################################################
# 'CN' = carbon-nitrogen single bond

# 'CNC' = carbon-nitrgen-carbon via single bonds
# 'CNO' = carbon-nitrgen-oxygen via single bonds
# 'CNN' = carbon-nitrgen-nitrogen via single bonds
# 'CNS' = carbon-nitrgen-carbon via single bonds

# 'CN(C)C' = tertiary amine
# 'CN(C)O' = hydroxylamine (OH-N) where N is bonded to two carbons
# 'CN(C)N' = hydrazine (N-N) where N is bonded to two carbons
# 'CN(C)S' = sulfenamide (N-S) where N is bonded to two carbons
# 'NC=O' = amide

# 'CO' = carbon-oxygen single bond

# 'COC' = carbon-oxygen-carbon via single bonds
# 'COO' = carbon-oxygen-oxygen via single bonds
# 'CON' = carbon-oxygen-nitrogen via single bonds
# 'COS' = carbon-oxygen-carbon via single bonds

# 'O=CO' = carbonyl with carbon bonded to oxygen via single bon
# 'COC=O' = ester
# 'O=COO' = ester where single bonded oxygen is bonded with another oxygen
# 'NOC=O' = ester where single bonded oxygen is bonded with a nitrogen
# 'O=COS' = ester where single bonded oxygen is bonded with sulfur
################################################################

features = ['CN','CNC','CNO','CNN','CNS','CN(C)C','CN(C)O','CN(C)N','CN(C)S','NC=O',
            'CO','COC','COO','CON','COS','O=CO', 'COC=O','O=COO','NOC=O','O=COS'] 


# give a dictionary of features and all the sub-features that should not count
# toward structure matches--i.e., in the case of COC, any similar CO submatches
# should not be counted toward the COC total count.
# Rule: given a feature and its sub-feature dependencies, for all combinations
# of the set of molecules within those sub-feature dependencies, if they match
# the given feature then the greater feature should not be counted.
# feature_dependency = {"COC": ["CO"],
#                       }


def get_molecule_substructure_match_matrix(molecules, substructures):
    all_matches = pd.DataFrame()
    for molecule in molecules:
        print("\t" + str(molecule))
        for pattern in features:
            print("\t\t" + str(pattern))
            substruct_matches = submatches(molecule, pattern)
            for substruct_match in substruct_matches:
                new_match = pd.Series([molecule, pattern, substruct_match],
                                      index = ['Canonical_SMILES','pattern','substruct_match'])
                all_matches = all_matches.append(new_match, ignore_index=True)
    # print(all_matches.head())
    return all_matches



def matches_count_pivot(all_matches):
    all_matches = pd.pivot_table(all_matches,
                                 values = ['substruct_match'],
                                 index = ['Canonical_SMILES', 'pattern'],
                                 columns = ['substruct_match'],
                                 aggfunc = len).reset_index()
    return all_matches



def is_duplicate_struct_from_substructs(target_substruct, all_substruct_matches):
    # only look at the combinations within all_substruct_matches
    # that doesn't include the target_substruct
    all_substruct_matches.remove(target_substruct)
    
    # if all molecules are the same between the target_substruct
    # and the unique list of all molecules in the combination of 2 patterns
    # from all_substruct_matches, then this is a duplicate
    
    # Use set comparisons for the order of the atom locations to not matter, and use
    # list comparisons for the order of the atom locations to matter
    found_duplicate = False
    for substruct_match_1 in all_substruct_matches:
        for substruct_match_2 in all_substruct_matches:
            if substruct_match_1 != substruct_match_2:
                # wait = input("Press enter to continue...")
                substruct_atom_set = set(list(substruct_match_1).extend(list(substruct_match_2)))
                if substruct_atom_set == set(target_substruct):
                    found_duplicate = True
                else:
                    found_duplicate = False
    # if there are any found duplicates, found_duplicate will be True, else False
    return found_duplicate


def remove_duplicate_matches(all_matches):
    all_matches_filtered = pd.DataFrame()
    for molecule_smiles, substruct_matches in all_matches.groupby(by=['Canonical_SMILES']):
        # for each pattern, get all the patterns that could be inflating its
        # score by looking it up in feature_dependency. If there are
        # any sets out of all the combinations of substruct matches that
        # are equal to one of the substruct matches in this pattern, do not
        # append to the final filtered all_matches_filtered (i.e., filter it out
        # by leaving it out)
        molecule_patterns_matched = substruct_matches['pattern'].unique()
        for pattern in molecule_patterns_matched:
            this_patterns_substruct_matches = list(substruct_matches[ substruct_matches['pattern'] == pattern]['substruct_match'].unique())
            for target_substruct in this_patterns_substruct_matches:
                duplicate_pattern = is_duplicate_struct_from_substructs(target_substruct = target_substruct,
                                                                        all_substruct_matches = this_patterns_substruct_matches)
                if not duplicate_pattern:
                    keep_pattern_match = pd.Series([molecule_smiles, pattern, target_substruct],
                                                   index = ['Canonical_SMILES','pattern','substruct_match'])
                    all_matches_filtered = all_matches_filtered.append(keep_pattern_match,
                                                                       ignore_index = True)
    return all_matches_filtered


def process_molecule(molecule, substructures):
    # Capture this whole function in a try/except block. If one molecule out of the
    # 999,999 we're processing throws an error, which could happen for a number of
    # reasons, we don't want the whole program to crash; just handle the error
    # and skip that molecule.
    try:
        # Print out the molecule's index into the terminal so we can see the progress
        # of the parallel processing :) you will see each index print out of order,
        # but at least it will give a general idea of where the processing is at;
        # about how many molecules it has processed. Note, you can make a global
        # variable to represent the total number of molecules processed, but this
        # would require some work figuring out how all the multiple processes running
        # at once would write to the global variable with out locking it (probably
        # send to a queue with a single writer), but this is good enough!
        molecule_ix = molecule_ix_dict[molecule]
        print(str(molecule_ix))
        
        # process all the possible pattern matches for this molecule
        all_matches = pd.DataFrame()
        for pattern in features:
            # print("\t\t" + str(pattern))
            substruct_matches = submatches(molecule, pattern)
            for substruct_match in substruct_matches:
                new_match = pd.Series([molecule, pattern, substruct_match],
                                      index = ['Canonical_SMILES','pattern','substruct_match'])
                all_matches = all_matches.append(new_match, ignore_index=True)
        # print(all_matches.head())
        
        # Save the molecule's pattern matches to CSV. Don't name the CSV the name of the
        # molecule because there could be special characters in the molecule's SMILES
        # which will throw a file system error.
        # Later we can run a script to read in all the CSVs into one giant DataFrame
        # or push to a database if you'd like to have the data that way.
        # Otherwise, you could write a script to iterate through each of the molecule's
        # CSV files to do the filtering of redundant pattern matches that way.
        # Probably have these filtered files in a different folder so you don't overwrite
        # this data, which will allow you to do trial & error developing the filtering routine.
        all_matches.to_csv(f"output/molecule_pattern_matches/{molecule_ix}.csv",
                           index = False)
    except Exception as e:
        # Log any errors, and also record the SMILES string which threw the error.
        # 
        log.error(f"{molecule}: {e}")
        


def main():
    # Process the molecule pattern matches in parallel to speed up execution time.
    # As soon as this runs you should see 1 CSV file pop up in the folder
    # output/molecule_pattern_matches for each molecule processed!
    pool = mp.Pool(mp.cpu_count())
    pool.map(partial(process_molecule, # function to process in parallel using multiple cores on your machine at once
                     substructures=features), # keyword argument(s) (extra arguments the function expects) here
             molecules # iterable to pass to function, where each entry is passed to the first argument of the function
             )
    
    # process serially if fast enough:
    # print("\tGenerating all possible combinations and saving the data...")
    # all_matches = get_molecule_substructure_match_matrix(molecules = chem_df['Canonical_SMILES'].unique(),
    #                                                      substructures = features)
    # all_matches.to_csv("all_matches.csv", index=False)
    # 
    # # pivot for easy visualizing of similar substruct matches across patterns
    # print("\tPivoting the data by pattern match result and saving...")
    # all_matches_pivot = matches_count_pivot(all_matches)
    # all_matches_pivot.to_csv("matches_pivoted.csv", index=False)
    # 
    # # remove pattern matches where the count is inflated by a smaller structure's
    # # pattern match
    # print("\tRemoving duplicates on the pattern matches and saving...")
    # filtered_matches = remove_duplicate_matches(all_matches)
    # filtered_matches.to_csv("filtered_matches.csv", index=False)

if __name__ == "__main__":
    main()

# FUNCTION FOR FINDING SUBSTRUCTURES IN MOLECULES

# For each of the functional groups in 'features' list, create a column for that functional group that counts how many times
# they occur in a given molecule. 

# for n in features:
#     df[n]=[len(submatches(i, n)) for i in df['Chem']]


# COUNTING ALCOHOL GROUPS (COH)

# The 'CO' column counts ethers/carboxylic acids as well as alcohols. 
# To remove ether contributation to this total, I can substract counts. There are 2 COs in 1 COC, hence the why we are 
# doubling df['COC'] below. 
# 
# df['COH'] = df['CO'] - (2*df['COC'] + df['O=CO']) 
# df.head()
# 
# 
# Chem.MolFromSmiles('COC')
# 
# 
# # COUNTING CARBOXYLIC ACID GROUPS (COOH)
# 
# # The 'O=CO' column may also count esters. Have to subtract these ester types from the count. 
# df['COOH'] = df['O=CO'] - (df['COC=O'] + df['O=COO'] + df['NOC=O'] + df['O=COS']) 
# df.head()
# 
# 
# 
# # COUNTING PRIMARY AMINE GROUPS (NH2)
# 
# # The substructure "NC" include any single carbon-nitrogen bond, including those
# # in secondary and teritary amines. If a molecule contains at least 1 teritary amine 'CN(C)C' as in Condition 1, then 
# # 3 NCs are counted. If any secondary amines 'CNC' exist, then 2 NCs are added to the count. The idea behind Condition 1 is
# # to subtract these "false alarms," leaving only the primary amines. Condtion 2 is the case where no tertiary amines
# # are in the molecule. 
# 
# df['NH2'] = '' # Create an empty column for primary amine groups
# df.loc[df['CN(C)C'] > 0, 'NH2'] = df['CN'] - (3*df['CN(C)C'] + 2*df['CNC'] + df['NC=O']) # Condition 1 for primary amine calculation
# df.loc[df['CN(C)C'] == 0, 'NH2'] = df['CN'] - (2*df['CNC'] + df['NC=O'])  # Condtion 2 for primary amine calculation 
# df.head()
# 
# 
# 
# # COUNTING SECONDARY AMINE GROUPS (C-NH-C)
# 
# 
# 
# # Polymerization Booleans
# df.loc[df['COH']>0 & df['COOH']>0, 'amide'] = 1
# df.loc[df['COH']==0 or df['COOH']==0, 'amide'] = 0



