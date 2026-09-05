# glucose, third generation

The glucose degradation network truncated at generation 3.

- `Data_cleaning.ipynb` — reads the network, strips the MØD rule syntax from the reaction names, orders species by first occurrence, pickles the two sides out
- `reactants`, `products` — Python pickles, each a dict from `reaction_0`, `reaction_1`, … to the SMILES on that side; together the stoichiometric matrix
- `Ematrix.m` — finds the extreme currents, the positive solutions of `SC·E = 0`. Guy Schmitz's 2008 program, with Kolar-Anic, Anic and Cupic; the same copy is in each folder
- `test.mat` — the assembled matrices as MATLAB reads them
