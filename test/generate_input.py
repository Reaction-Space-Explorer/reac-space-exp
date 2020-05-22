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


def generate_new_input(current_file, new_file, reaction_model, include_edge=False):
    new_input = open(new_file, 'w')
    # remove simpleReactor from prev file
    raw_code = clean_prev(current_file)
    new_code = []
    inerts = ['He', 'Ne', "Ar"]
    for sp in reaction_model.core.species:
        # The inert gases which RMG adds by default must not be saved in the input file
        if sp.label not in inerts:
            new_code.append('\n')
            new_code.append('species(')
            new_code.append("\tlabel='{0}',".format(sp.label))
            new_code.append('\treactive=True,')
            new_code.append(f"\tstructure=SMILES('{sp.smiles}')")
            new_code.append(')')
    # if include_edge:
    # do the same for edge species
    # now add a simpleReactor()

    new_code.append('\n')
    new_code.append('simpleReactor(')
    new_code.append(f"\ttemperature=({temperature},'K'),")
    new_code.append(f"\tpressure=({pressure},'bar'),")
    new_code.append(f"\tterminationTime=(1e0,'s'),")  # by default make all reactions run for 1s
    new_code.append("\tinitialMoleFractions={")

    # I'm assuming we won't add edge species
    # We'll have to expand this method if we do so

    # let all species have equal mole fractions initially
    # # i did a -3 there because core.species contains He, Ne and Ar which we wanna exclude
    average_mol_frac = 1.0 / (len(reaction_model.core.species)-3)
    for sp in reaction_model.core.species:
        if sp.label not in inerts:
            new_code.append(f"\t\t'{sp.label}' : {average_mol_frac},")
            if reaction_model.core.species.index(sp) == (len(reaction_model.core.species) - 1):
                new_code[len(new_code) - 1].replace(',', '')  # remove the comma for the last item
    new_code.append('\t}')
    new_code.append(')')

    # finally, write all of this
    new_input.writelines(raw_code)  # lines from the old input_0.py don't need an \n added
    for line in new_code:
        new_input.write('\n')
        new_input.write(line)

    new_input.close()