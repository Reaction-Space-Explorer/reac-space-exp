# radical network

The network from [`radicals/`](../../radicals), the same treatment as the glucose networks and far smaller.

- `Data_cleaning.ipynb` — reads the network, strips the MØD rule syntax from the reaction names, orders species by first occurrence, pickles the two sides out
- `reactants`, `products` — Python pickles, each a dict from `reaction_0`, `reaction_1`, … to the SMILES on that side; together the stoichiometric matrix
- `Ematrix.m` — finds the extreme currents, the positive solutions of `SC·E = 0`. Guy Schmitz's 2008 program, with Kolar-Anić, Anić and Čupić (J. Phys. Chem. A 112, 13452, doi:10.1021/jp8056674); the same copy is in each folder
- `test.mat` — the assembled matrices as MATLAB reads them
