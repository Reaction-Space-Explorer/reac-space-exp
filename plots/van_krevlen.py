import matplotlib.pyplot as plt

"""
This is a MOD script to create the van Krevlen plot
run this using ```mod -f van_krevlen.py```
"""

# store the generation in which the molecule appeared
gens_list = []
# list within a list; each item within this will hold a list containing 3 elements: number of C, H and O.
atom_count_list = []
# count C, H and O
atoms_of_interest = "CHO"
# list of computed ratios
hc_ratio_list = []
oc_ratio_list = []

# open the output .txt
with open('../main/glucose/glucose_degradation_output.txt') as glu_out:
    contents = glu_out.readlines()
    for line in contents:
        # get rid of \n at the end of the line
        line = line.rstrip('\n')
        # break line into components separated by a tabspace
        comps = line.split("\t")
        gen_num = int(comps[0][1])
        mol = smiles(comps[1])
        atom_count = []
        for atom in atoms_of_interest:
            atom_count.append(mol.vLabelCount(atom))
        gens_list.append(gen_num)
        atom_count_list.append(atom_count)
        # also calculate the ratio of atom numbers right here?
        hc_ratio = atom_count[1]/atom_count[0]
        oc_ratio = atom_count[2]/atom_count[0]
        hc_ratio_list.append(hc_ratio)
        oc_ratio_list.append(oc_ratio)


# make the figure container
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

# and color by generation
colors_for_gen = ['orange', 'steelblue', 'greenyellow', 'pink', 'violet']

# now, for each generation
for gen in range(max(gens_list), 0, -1):
    hc_ratios = [hc_ratio_list[i] for i in range(len(hc_ratio_list)) if gen == gens_list[i]]
    oc_ratios = [oc_ratio_list[i] for i in range(len(oc_ratio_list)) if gen == gens_list[i]]
    plt.plot(oc_ratios, hc_ratios, color=f'{colors_for_gen[gen-1]}', linestyle='', marker='o',
                alpha=0.25, label=f'Generation {gen}')
plt.legend(loc='upper left')
plt.xlabel('O:C Ratio')
plt.ylabel('H:C Ratio')
plt.title('van Krevlen diagram for the glucose network')
#fig.tight_layout()
plt.savefig('van_krevlen_diag.jpg', dpi=150)
plt.show()