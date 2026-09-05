# Stoichiometric network analysis

The generated networks read as stoichiometric matrices, by **Alejandro Lozano**
([@alejandro-lozano-dev](https://github.com/alejandro-lozano-dev)), imported here with its
history.

- `glucose_3rd_generation/` — glucose degradation at generation 3
- `glucose_5th_generation/` — at generation 5, with species matched and missed against the reference list
- `radical_network/` — the network from [`radicals/`](../radicals)
- `glucose_degradation_output.txt` — the network the analyses read
- `lectures/`, `course_homework/` — third-party course material on the method, not analyses

Each folder holds the same four files: `Data_cleaning.ipynb` builds the matrix and pickles it
into `reactants` and `products`, `Ematrix.m` finds the extreme currents, `test.mat` holds the
result. `Ematrix.m` is Guy Schmitz's 2008 program, not written here.

`lectures/` and `course_homework/` came with the original repository and are third-party
teaching material rather than work on these networks: the slides are Bernhard Palsson's, from
*Dynamic States of Biological Networks*, and the homework is that course's. Each folder says
so. [MassPy](https://masspy.readthedocs.io/en/latest/notebooks/SB2_textbook/SB2-Chapter-1.html)
publishes the textbook notebooks openly.
