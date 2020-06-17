from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# molVs' smiles might not be identical to rdkit's so let's use RDKit instead of molVs
enum = rdMolStandardize.TautomerEnumerator()


def clean_taut(dg, dg_execute):
    subset = dg_execute.subset
    universe = dg_execute.universe
    # The following line may seem redundant but it is needed for retaining the indices of subset
    graphs = [g for g in subset]
    # list  of smiles associated with each Graph object
    mod_subset_smiles = [g.smiles for g in subset]
    # a mapping of the Graph objects with their corresponding DGVertex objects
    mod_dgverts = {v.graph: v for v in dg.vertices if v.graph in subset}
    # the list of RDKit canonical SMILES
    rdkit_subset_smiles = []

    taut_tautset_dict = {}
    # avoid duplicates in to_remove
    to_remove = set()
    # Populate all possible tautomers and add to the above empty dict
    for mod_smiles in mod_subset_smiles:
        mol = Chem.MolFromSmiles(mod_smiles)
        # convert into rdkit canonical smiles
        rdkit_smiles = Chem.MolToSmiles(mol)
        rdkit_subset_smiles.append(rdkit_smiles)
        # generate Mol objects of all possible tautomers
        all_tauts = enum.Enumerate(mol)
        # convert into smiles
        all_taut_smiles = tuple(Chem.MolToSmiles(taut) for taut in all_tauts)
        #print(all_taut_smiles)
        taut_tautset_dict[mod_smiles] = all_taut_smiles

    list_tautset = tuple(taut_tautset_dict.values())
    class_ids = {}
    # create initial class ids
    for i in range(len(list_tautset)):
        class_ids.update({list_tautset[i]: i})

    for i in range(len(list_tautset)):
        for j in range(i+1, len(list_tautset)):
            for item in list_tautset[i]:
                if item in list_tautset[j]:
                    #print(f'Classes {i} and {j} have common elements')
                    #class_ids[list_tautset[i]] = class_ids[list_tautset[i]]
                    class_ids[list_tautset[j]] = class_ids[list_tautset[i]]
                    #print('Updated class ids to ', class_ids[list_tautset[i]])
                    break
    # Now, since we're done finding identical class ids
    # A dictionary of observed tautomer classes mapped as {class id: [list of smiles]}
    dict_observed_tauts = {}
    for taut_class, class_id in class_ids.items():
        for i in range(len(rdkit_subset_smiles)):
            if rdkit_subset_smiles[i] in taut_class:
                # the smiles of the mod graph
                if class_id in dict_observed_tauts.keys():
                    dict_observed_tauts[class_id].append(mod_subset_smiles[i])
                else :
                    dict_observed_tauts[class_id] = [mod_subset_smiles[i]]
    print('{0} unique tautomer classes in {1} molecules'.format(len(dict_observed_tauts.keys())
                                , len(subset)))
    for taut_class in dict_observed_tauts.values():
        if len(taut_class) > 1:
            for i in range(1, len(taut_class)):
                to_remove.add(taut_class[i])
    for smiles in to_remove:
        # The DGVertex associated with the molecule to be removed
        dg_vertex = mod_dgverts[graphs[mod_subset_smiles.index(smiles)]]
        # TODO: d.right should point towards the graph that was preserved
        for e in dg_vertex.inEdges:
            for source in e.sources:
                    d = Derivations()
                    d.left = [source.graph]
                    d.rules = e.rules
                    d.right = [graphs[mod_subset_smiles.index(smiles)]]
                    b.addDerivation(d)
        subset.remove(graphs[mod_subset_smiles.index(smiles)])
        universe.remove(graphs[mod_subset_smiles.index(smiles)])
    return subset, universe