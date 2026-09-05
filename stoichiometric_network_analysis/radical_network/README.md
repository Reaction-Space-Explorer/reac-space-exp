# radical network

The radical network from [`radicals/`](../../radicals) read as a stoichiometric matrix, the
same treatment as the glucose networks and far smaller.

| | |
|---|---|
| `Data_cleaning.ipynb` | Reads the network, strips the MØD rule syntax from the reaction names, orders the species by first occurrence and pickles the two sides out |
| `reactants`, `products` | Python pickles, each a dict from `reaction_0`, `reaction_1`, … to the SMILES on that side of it. Together they are the stoichiometric matrix |
| `Ematrix.m` | Finds the extreme currents, the positive solutions of `SC·E = 0`. Not written here: it is Guy Schmitz's 2008 program, with Ljiljana Kolar-Anic, Slobodan Anic and Zeljko Cupic, from *Stoichiometric network analysis and associated dimensionless kinetic equations*. The same copy sits in each folder |
| `test.mat` | The assembled matrices, as MATLAB reads them |
