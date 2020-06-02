"""
Contains methods for generating new input_0.py files
This will help us do more "generations" of reactions
"""

# I couldn't find a way to access these from the reaction model itself
# TODO: Access these from the models, instead of providing it explicitly
temperature = 800
pressure = 1.0


def clean_prev(current_file):
    """
    Intended to "clean" (remove) simpleReactor() from an input_0.py
    :param current_file: the file that was just used for a run
    :return: code without simpleReactor()
    """
    prev_input = open(current_file)
    prev_code = prev_input.readlines()
    prev_input.close()
    # I will be removing simpleReactor() lines that contained in the previous file
    starting_line = 0
    end = 9999
    found = [False, False]
    while not found[1]:
        for line in prev_code:
            if "simpleReactor" in line:
                found[0] = True
                starting_line = prev_code.index(line)
            if found[0] and line == ')':
                found[1] = True
                end = prev_code.index(line)
                break
        if not found[1]:
            raise ValueError('Cannot find simpleReactor(), please check input file')
    del prev_code[starting_line:end + 1]  # the last line needs to be included hence end+1
    return prev_code


def generate_new_input(current_file, new_file, reaction_model, old_species_labels,
                       include_edge=False):
    new_input = open(new_file, 'w')
    # remove simpleReactor from prev file
    raw_code = clean_prev(current_file)
    new_code = []
    inerts = ['He', 'Ne', "Ar"]
    species = reaction_model.core.species
    if include_edge:
        species = species + reaction_model.edge.species
    # remove list of species in previous generation
    # because they already exist in the input file
    species = [sp for sp in species if sp.label not in old_species_labels]
    species = [sp for sp in species if sp.label not in inerts]
    '''for sp in species:
        # The inert gases which RMG adds by default must not be saved in the input file
        if sp.label in old_species_labels or sp.label in inerts:
            print('Removing ' + sp.label)
            print('List had ' + str(len(species)) + ' species before')
            species.remove(sp)
            print(str(len(species)) + ' now')'''
    for sp in species:
        new_code.append('\n')
        new_code.append('species(')
        # Note: RMG doesn't allow labels to have '+' sign
        # So don't treat labels as exact SMILES strings
        new_code.append("\tlabel='{0}',".format(sp.label.replace('+','')))
        new_code.append('\treactive=True,')
        new_code.append(f"\tstructure=SMILES('{sp.smiles}')")
        new_code.append(')')
    # now add a simpleReactor()

    new_code.append('\n')
    new_code.append('simpleReactor(')
    new_code.append(f"\ttemperature=({temperature},'K'),")
    new_code.append(f"\tpressure=({pressure},'bar'),")
    new_code.append(f"\tterminationTime=(1e0,'s'),")  # by default make all reactions run for 1s
    new_code.append("\tinitialMoleFractions={")

    # let all species have equal mole fractions initially
    average_mol_frac = 1.0 / (len(species))
    for sp in species:
        if sp.label not in inerts:
            new_code.append(f"\t\t'{sp.label.replace('+','')}' : {average_mol_frac},")
            if species.index(sp) == (len(species) - 1):
                new_code[len(new_code) - 1].replace(',', '')  # remove the comma for the last item
    new_code.append('\t}')
    new_code.append(')')

    # finally, write all of this
    new_input.writelines(raw_code)  # lines from the old input_0.py don't need an \n added
    for line in new_code:
        new_input.write('\n')
        new_input.write(line)
    new_input.close()