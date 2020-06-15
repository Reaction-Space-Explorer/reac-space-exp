# Installation
I made this tut because I thought some people were having trouble setting MOD up. I made this guide for Ubuntu based builds but I think you can build it on Windows too. A Ubuntu VM might work equally well.

Some of the steps might seem like a copy-off from the installation website but I tried to add a few words to make life easy for others. I've also updated to provide resolutions for many bugs.

First off, install the following (which is a huge list)
* CMake
    * You can install the GUI version if you'd like (it makes life easier). Note that on Ubuntu (at least on my version), the GUI doesn't run for some reason but even when it doesn't, one can make full use of the command-line way of doing things without any trouble.
    * The following should install it in Linux
    ```bash
    sudo apt install make
    ```
* A C++ Compiler with C++14 support. As per the installation website, GCC version 6.1+ should work.
    ```bash
    sudo apt install g++ gcc
    ```
* Boost with Boost.
    * **UPDATE:** People with Ubuntu 20.04 should install Boost using **conda**.
    * Some people in general had trouble installing using source code. Try using the following commands instead
    ```bash
    sudo apt install -y libboost-dev libboost-python-dev libboost-all-dev
    ```
    * To build using source code: Download the .tar.gz containing the source from [here](https://www.boost.org/)
    * After unpacking, open it and use the following bash commands
    ```bash
    ./bootstrap.sh --with-python=python3
    ./b2
    ./b2 install
    ```
* Graphviz
    * Update: the distribution for Ubuntu doesn't come with RSVG or Pangocairo. Install them using the commands
    ```bash
    sudo apt install -y librsvg2-dev libpango1.0-dev pdf2svg
    ```
    * After you've installed the above two, you can now build Graphviz from the [source](https://graphviz.gitlab.io/_pages/Download/Download_source.html). The link has listed the terminal commands you'll need as well.
    * Note that you need CMake installed for building this.
* Some LaTeX distribution (MOD uses LaTeX to generate the PDF output) with some science packages (they were needed in my case)
    * Here are some terminal commands to make it easier for you
    ```bash
    sudo apt install -y texlive-base tex-livescience texlive-latex-base texlive-latex-extra
    ```
* OpenBabel
    ```bash
    sudo apt-get install -y openbabel libopenbabel-dev
    ```
    * Note: Apparently MOD wants to use some version of 2.0 series. If "babel -v" prints something like OpenBabel 2.3.2, there should be no problem. There was a problem for those who had 3.0.
    * Note 2: In Ubuntu 20.04, older versions of most packages aren't available under ```sudo apt``` (so this applies even to openbabel and libopenbabel-dev as well) but this can be resolved. There is a file named **sources.list** in the folder ```/etc/apt/``` which contains the list of sources where ```apt``` gets its list of packages from. Add this line to the bottom of the file and save
    ```bash
    deb http://cz.archive.ubuntu.com/ubuntu eoan main universe
    ```
    Then run the following in the terminal
    Now you should be able to access **eoan** packages(Eoan Ermine is the name of Ubuntu 19.10 release), which supported the versions we desire. Don't worry, everything's safe.

    Then you need to do ```sudo apt update``` which will refresh the list of packages available (this time, it will include those eoan packages). You can then do:
    ```bash
    sudo apt install openbabel=2.4.1+dfsg-3 libopenbabel-dev=2.4.1+dfsg-3
    ```
    Note that we have told **apt** to forcefully install specific versions of these packages.
* Sphinx (which is a documentation generator for Python)
    * I think I used the command-line
    ```bash
    sudo apt install -y python3-sphinx
    ```

## Clone the GitHub repository
Once you have all those dependencies installed (which would take a bit of time but there's no alternative), clone the [GitHub repo](https://github.com/jakobandersen/mod) of MOD

```bash
git clone https://github.com/jakobandersen/mod.git
```
## Some configs and then make build directories
Important: this will fetch some dependencies like graph_cannon (created by Jakob Andersen himself). To be able to do this, you need to have cloned using Git. Open the folder in which you cloned MOD repo 
```bash
git submodule update --init --recursive
./bootstrap.sh
```

After this, create a directory named build and make the build files usinsg CMake
```bash
mkdir build
cd build
cmake ../ <options>
make -j 
make install
```
In the ```<options>```, you **have to** turn -DBUILD_PY_MOD=on and -DBUILD_POST_MOD=on for our purposes. I recommend passing almost all options (atleast those which were =on) on the main [installation page](http://jakobandersen.github.io/mod/installation.html). For your convenience, here's what I used

```bash
cmake ../ -DBUILD_DOC=on -DBUILD_POST_MOD=on -DBUILD_PY_MOD=on -DBUILD_TESTING_SANITIZERS=on -DENABLE_SYMBOL_HIDING=on -DENABLE_DEP_SYMBOL_HIDING=on -DENABLE_IPO=on -DUSE_NESTED_GRAPH_CANON=on -DWITH_OPENBABEL=on
```

*Note: The above process takes quite a bit of time.*
In case you get some error like "Can't write/create" or "Permission denied" or anything of the sort in the last command. Try using **sudo** prefix, which is intended to give you root privileges while executing the command

```bash
sudo make install
```
This fixed it for me
## How to tell if I installed MOD properly
Type the following command in the terminal
```bash
mod --version
```
This should give you an output telling you the prefix, version, etc. of the current installation of mod. If this was installed properly, we're all set.
### **Note**:
Some people on Ubuntu 20.04 (including myself) or those with problems with Boost received a weird error while trying to run **mod**.
```python
ImportError: /usr/lib/x86_64-linux-gnu/libboost_python38.so.1.67.0: undefined symbol: _Py_tracemalloc_config
```
Apparently ```sudo``` prefix before the **mod** command seemed to have had fixed it. But that's an illusion: the problem is that the boost libraries were built with python 3.8 (Ubuntu 20.04 comes with 3.8) and a different python version is being used to run it. This is because ```sudo``` overrides environment variables and aliases (check the difference between ```python3 --version``` and ```sudo python3 --version```)
Don't use ```sudo``` to run any code (as Dr. Andersen himself said). The solution for me was to install Boost using **conda**.

```bash
conda install -c anaconda boost
```
And pass an additional option to **cmake** to force it to use this version of boost ```-DBOOST_ROOT=~/anaconda3```

## How to run files with MOD
If you have a file you want to run, type the command
```bash
mod -f test.py 
```

Now, you might want to be able to import the mod package by writing
```python
import mod
```
in your code. But as you might figure out that importing it gives you a package not found error. As the solution, you need to add the folder containing the Python bindings of MOD to your **PYTHONPATH**
. To do this, find the "prefix" by running 
```bash
mod
```
which outputs something like
```bash
MØD Wrapper --------------------------------------------------------
Prefix: /usr/local/bin/..
MedØlDatschgerl version 0.9.0.1
Plugins ------------------------------------------------------------
Global (/usr/local/bin/../share/mod/plugins):
 mod: /usr/local/lib
```
As you can see, for me, the prefix was **/usr/local/bin**. Now the Python bindings aren't in the /bin/ folder but in the /lib/ folder.
So what you want to add to your PYTHONPATH is **/usr/local/lib/**
One problem is that if you add to your PYTHONPATH using a terminal in Ubuntu, the addition will be reset once you close that terminal.

For a more permanent addition, go to your $HOME directory and view the hidden files. You'll find a file named **.bashrc** which is a script that is run everytime a terminal is created. Add this line to the bottom of the file.
```bash
export PYTHONPATH=$PYTHONPATH:/usr/local/lib/
```
## Still having trouble?
Some of your problems might be well known, so google your errors and see if you can fix them. Message me on Slack in the #reactionnetworkgeneration channel or personally (whichever you prefer).