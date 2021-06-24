# Installation
I made this tut because I thought some people were having trouble setting MOD up. I made this guide for Ubuntu based builds but I've heard it can be built on [Windows subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) as well. A Ubuntu VM on Windows might work equally well. And personally, I now use it on Manjaro which is based on Arch.

**2021 Update**: A ```conda``` package as well as a Docker image is available for MOD (see [this page](http://jakobandersen.github.io/mod/installation.html)), though for a few people, conda couldn't do the job due to some weird reasons. A lot of the dependencies I listed below can be installed quickly if you do it using the ```bindep``` file (see [official compilation from source instructions](http://jakobandersen.github.io/mod/compiling.html)) and the also the ```requirements.txt```. Nevertheless, I believe the following is still relevant in case someone gets stuck. Let me know if something is outdated and I will try to fix that.

These are based on the instructions on Andersen's website but I tried to add a few words and resolution of a few errors to make life simpler for others.

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
    sudo apt install -y texlive-base texlive-science texlive-latex-base texlive-latex-extra
    ```
* OpenBabel
    ```bash
    sudo apt-get install -y openbabel libopenbabel-dev
    ```
    * *2021 Update*: MOD works fine with openbabel v3
#### Note 2:
In Ubuntu 20.04, certain versions of a few packages aren't available under ```apt``` (so this applies even to openbabel and libopenbabel-dev as well) but this can be resolved. There is a file named **sources.list** in the folder ```/etc/apt/``` which contains the list of sources where ```apt``` gets its list of packages from. To edit this file, you need root privileges, so perhaps running the text editor using ```sudo``` should help.
    ```sudo gedit /etc/apt/sources.list```
    Then add this line to the bottom of the file and save
    ```bash
    deb http://cz.archive.ubuntu.com/ubuntu eoan main universe
    ```
    Then run the following in the terminal
    Now you should be able to access **eoan** packages (Eoan is the name of Ubuntu 19.10 release), which supported the versions we desire. Don't worry, everything's safe.

    Then you need to do ```sudo apt update``` which will refresh the list of packages available (this time, it will include those eoan packages). You can then do:
    ```bash
    sudo apt install openbabel libopenbabel-dev
    ```
    Note that we have told **apt** to forcefully install specific versions of these packages.
* **Optional:**Sphinx, which is a documentation generator for Python, is not needed if you're setting ```-DBUILD_DOC=off``` while running ```cmake``` to build MOD. Personally, I never used the local documentation anyways and relied on 
    * I think I used the command-line
    ```bash
    pip install -U sphinx
    ```
    * Again, ```apt``` might not have the latest version for sphinx, I think ```pip``` is safer.


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
In the ```<options>```, you **have to** turn ```-DBUILD_PY_MOD=on``` and ```-DBUILD_POST_MOD=on``` for our purposes. Check the main [compilation from source instructions](http://jakobandersen.github.io/mod/compiling.html) for detail about these . For your convenience, here's what I used

```bash
cmake ../ -DBUILD_DOC=off -DBUILD_POST_MOD=on -DBUILD_PY_MOD=on -DBUILD_TESTING_SANITIZERS=on -DENABLE_SYMBOL_HIDING=on -DENABLE_DEP_SYMBOL_HIDING=on -DENABLE_IPO=on -DUSE_NESTED_GRAPH_CANON=on -DWITH_OPENBABEL=on
```
The log of this might be worth keeping for inspection, in case you want to see which version of which dependency it picked up. The next step is to
```
make
```
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
Apparently ```sudo``` prefix before the **mod** command seemed to have had fixed it. But that's an illusion: the problem is that the boost libraries were built with python 3.8 (Ubuntu 20.04 comes with 3.8) and a different python version is being used to run it. This is because ```sudo``` overrides environment variables and aliases (check the difference between ```python3 --version``` and ```sudo python3 --version```). The solution for me was to install Boost using **conda**.

```bash
conda install -c anaconda boost
```
And pass this additional parameter to ```cmake``` to force it to search this installation of boost: ```-DBOOST_ROOT=~/anaconda3``` (it would be ```miniconda3``` if you are using that instead).

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
which outputs something like (note: this is a year old and the latest mod version is **0.12** but the essence is the same).
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
Some of your problems might be well known, so google your errors and see if you can fix them. Message me on Slack or personally (whichever you prefer).