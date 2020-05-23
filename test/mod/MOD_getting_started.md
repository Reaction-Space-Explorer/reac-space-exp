# Installations
I made this tut because I thought some people were having trouble setting MOD up. I made this guide for Ubuntu based builds but I think you can build it on Windows too. A Ubuntu VM might make life easy as well.

Some of the steps might seem like a copy-off from the installation website but I wrote them (with a few added words) to make life easy for others.

First off, install the following (which is a huge list)
* CMake
    * You can install the GUI version if you'd like (it makes life easier). Note that on Ubuntu (at least on my version), the GUI doesn't run for some reason but even when it doesn't, one can make full use of the command-line way of doing things without any trouble.
* A C++ Compiler with C++14 support. As per the installation website, GCC version 6.1+ should work.
* Boost with Boost.Python
    * Download the .tar.gz containing the source from [here](https://www.boost.org/)
    * After unpacking, open it and use the following bash commands
    ```bash
  ./bootstrap.sh --with-python=python3
  ./b2
  ./b2 install
    ```
* Graphviz
    * I think I had to build this from the [source](https://graphviz.gitlab.io/_pages/Download/Download_source.html). The link has listed the terminal commands you'll need as well.
    * Note that you need CMake installed for this.
* Some LaTeX distribution (MOD uses LaTeX to generate the PDF output) with some science packages (they were needed in my case)
    * Here are some terminal commands to make it easier for you (I'm unsure if the last one was needed for me)
```bash
sudo apt install texlive-core
sudo apt install texlive-science
sudo apt install texlive-latex-extra
```
* OpenBabel
```bash
sudo apt-get install libopenbabel-dev
```
* Sphinx (for the sake of documentation, I think)
    * I think I was able to install it using sudo apt-get

## Clone the GitHub repository
Once you have all those dependencies installed (which would take a bit of time but there's no alternative), clone the [GitHub repo](https://github.com/jakobandersen/mod) of MOD

```bash
git clone https://github.com/jakobandersen/mod.git
```
## Some configs and then make build directories
Open the folder in which you cloned MOD repo 
```bash
git submodule update --init --recursive
./bootstrap.sh
```

After this, create a directory named build and make the build files usinsg CMake
```nashorn js
mkdir build
cd build
cmake ../ <options>
make -j 
make install
```
In the options, you **have to** turn -DBUID_PY_MOD=on and -DBUILD_POST_MOD=on for our purposes. I recommend applying almost all (atleast those which were =on) on the main [installation page](http://jakobandersen.github.io/mod/installation.html)

*Note: The above process takes quite a bit of time.*
In case you get some error like "Can't write" or "Can't create" or anything of the sort in the last comman. Try using **sudo** prefix

```bash
sudo make install
```
This fixed it for me
## How to tell if I installed MOD properly
Type the following command in the terminal
```nashorn js
mod
```
This should give you an output telling you the prefix, version, etc. of the current installation of mod. If this was installed properly, we're all set.

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
You may use the following command in the terminal
```bash
export PYTHONPATH=$PYTHONPATH:/usr/local/lib/
```

## Still having trouble?
Message me on Slack in the #reactionnetworkgeneration channel or personally (whichever you prefer) but public has the advantage that someone else might be struggling with the same problem.