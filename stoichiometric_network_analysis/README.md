# Stoichiometric network analysis

The generated networks read as stoichiometric matrices, by **Alejandro Lozano**
([@alejandro-lozano-dev](https://github.com/alejandro-lozano-dev)), imported here with its
history.

- `glucose_3rd_generation/` — glucose degradation at generation 3
- `glucose_5th_generation/` — at generation 5, with species matched and missed against the reference list
- `radical_network/` — the network from [`radicals/`](../radicals)
- `glucose_degradation_output.txt` — the network the analyses read

Each folder holds the same four files: `Data_cleaning.ipynb` builds the matrix and pickles it
into `reactants` and `products`, `Ematrix.m` finds the extreme currents, `test.mat` holds the
result. `Ematrix.m` is Guy Schmitz's 2008 program, not written here.

The course lectures and worked homework that came with the original repository are not here,
being third-party teaching material. [MassPy](https://masspy.readthedocs.io/en/latest/notebooks/SB2_textbook/SB2-Chapter-1.html)
publishes the textbook notebooks openly.
