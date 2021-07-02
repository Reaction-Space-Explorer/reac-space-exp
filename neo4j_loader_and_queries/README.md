
# Welcome to the cycles pattern match sub-project!
![network pattern match example result](query_images/Size4_Pattern Match Rxn and Molecule Labels.JPG)


# Installation
Installation instructions below assume that you have already downloaded this
git repository, and you have the folder `reac-space-exp` on your drive somewhere,
ready to go!


## Python Installation
### Base Python
The script `main.py` was developed using [Python version 3.7.3 from python.org](https://www.python.org/downloads/release/python-373/).

For experienced users, see `requirements.txt` file for installation.

For more help, keep reading...

If you are installing Python it is recommended that you do not add Python to PATH, but
if you are sure that this is the only version of Python you are going to install
and you want to make the commands easier, then go ahead and check "Add to PATH".

Open a CMD prompt. If you are not familiar with CMD prompt, you type a command
and press Enter to run it. Be careful; there is no undo. Up/down arrows to access
previous commands you entered. Tab to auto-complete file names.

Get into this directory current directory by running the command:
```
> cd C:/path/to/wherever/you/have/folder/reac-space-exp/neo4j_loader_and_queries
```

Or, open File Explorer to this folder, and type `cmd` in the search bar--the CMD
prompt will open to the directory at this location.

If you added Python to PATH, this is how you run Python from CMD prompt:
```
> python
```
If you did not, search for the `python.exe` file. If you are using Windows 10,
I like to find it by hitting `win` key, typing `python`, right click on the Python
result and select `Open file location`, then right click on the Python shortcut
and again select `Open file location`. You should now be in the folder with `python.exe`.
Hold shift and right click on `python.exe`, select `Copy as path`. Now paste the
path in CMD prompt.
```
> "C:\Users\yourusername\AppData\Local\Programs\Python\Python37-32\python.exe"
```

Any time below I use "python" as a command in CMD prompt, you must use this full
path instead.

If you are running Python in the CMD prompt and can't figure out how to get out,
running this Python command will get you out of Python and back into CMD prompt:
```
>>> exit()
```
### Creating a Virtual Environment
If you plan on working on more than 1 Python project, it might be useful to create
a virtual environment to keep the projects separate. This step is not necessary
if not; you can skip to the next step Installing Packages.

Making sure you are still in the folder `neo4j_loader_and_queries` as your
current working directory, run:
```
> python -m venv venv
```
This is telling the virtual environment program to create a directory `venv`,
with a copy of your Python installation. Creating virtual environments by project
ensures that you don't have conflicting versions of libraries between projects,
since all the libraries are kept separate.

Your path to Python is now:
```
> "venv\Scripts\python.exe"
```

### Installing Packages
Using whichever Python you prefer (base installation or the virtual environment
Python if you've decided to do that), run the command:
```
> python -m pip install -r requirements.txt
```
This is telling Python's package manager pip to install all the libraries at
their specific version as specified in the `requirements.txt` file.


## Neo4j Installation
Download [Neo4j Desktop from neo4j.com](https://neo4j.com/download/?ref=try-neo4j-lp)
and install.



# Running the main.py Program
Run Neo4j. Click "Add Graph" to create a database. Give it a name, and password
`0000`. (If you'd like a different password, you can change it, just remember
to update `url = "bolt://neo4j:0000@localhost:7687"` in `main.py`).

On your new graph instance, click the Start button. Click Manage. Clicking Open
Browser will bring you to the Neo4j Browser, which you can use to do on-the-fly
CQL/cypher queries to the graph database and view the output. (This is how I
developed the main pattern match query in `graph_queries/_FINAL_QUERY_PARAMETERIZED.txt`.)

In `main.py`, make sure you have `EXPORT_PATHS` defined. An export path is a path
to a folder that includes data exported from MOD, which the folder `reac-space-exp/main`
in this project produces as the output. This data will detail the proliferation of
the chemical network at each generation, which `main.py` imports into Neo4j so
that it can run pattern match queries.

Using your path to Python that you installed the libraries on, do the following:
```
> python main.py
```

You should see the CMD prompt output something like this if it's working:
![example cmd prompt output](example_console_output.JPG)

If you get an error that Neo4j has run out of memory, stop the database instance,
then go to Manage > Settings. Try these memory settings (if your machine has enough
memory):
```
dbms.memory.heap.initial_size=18G
dbms.memory.heap.max_size=18G
dbms.memory.pagecache.size=18G
```
You can now restart the database instance and try running `main.py` again.

Each time the program runs, it will create a folder in the `output` directory,
and subfolders for each generation number it collected data on. Pattern match
results are stored in `autocat_query_results.csv`.


# Development tips
The script `main.py` starts at the bottom with the statement `if __name__ == "__main__":`.
I'd suggest to start reading there and following where it leads.

To develop on the pattern match query, view the files `graph_queries/example_query.txt`
and `_FINAL_QUERY_PARAMETERIZED.txt`. You can paste and run the full `example_query.txt`
in the Neo4j Browser (Manage > Open Browser), but for the parameterized version,
you have to replace `{{}}` variables with an actual value. The function
`get_tabulated_possible_autocatalytic_cycles()` reads in this parameterized text
file and replaces parameters with parameters given in the Python function. This 
makes it easy to programmatically change query parameters. `ring_size_range` and
`num_structures_limit` were set to limit how long it takes the query to run,
but feel free to modify these parameters.



# Output
The output folders are where the data from each run is collected, and the folder
is named by the timestamp of when it was created.

Below is an explanation of what each file in these folders means.

1. Neo4j_Imports Folder:
    This folder is copied from the source Neo4j_Imports folder and recorded here so that if you re-generate the reaction network with different parameters, you still have the network data from MOD from which this Neo4j network was generated.

2. query.txt:
    This file is the exact string passed to Neo4j to execute. The query parameters were passed programmatically, by the function call get_tabulated_possible_autocatalytic_cycles(). See also the file ../graph_queries/_FINAL_QUERY_PARAMETERIZED.txt. Query parameters are found and replaced using {{VARIABLE}} notation.

3. query_parameters.csv:
    This is a simple recording of the variables passed through get_tabulated_possible_autocatalytic_cycles() at the time the script was ran.

4. query_results.csv:
    This is the raw data output of the query, where each row represents a pattern matched by the query (where our target pattern is possible autocatalytic cycles). Has the same exact data as query_results.json, only it is more "human readable". To understand what the structure names in the data mean, refer to ../graph_queries/_query_visualization_diagram.jpg.

5. query_results.json:
    This is the raw data output of the query, where each index in the array represents a pattern matched by the query (where our target pattern is possible autocatalytic cycles). Has the same exact data as query_results.csv, only this format is easier to convert between dictionary/array objects. To understand what the structure names in the data mean, refer to ../graph_queries/_query_visualization_diagram.jpg.

6. ring_size_distribution.png:
    This plot is a simple plot of the number of ring patterns found by the number of molecules in the ring. Note that this is highly unreliable when matching a small number of patterns.



