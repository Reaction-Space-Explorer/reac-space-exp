# Using Ambit in Python
Ambit-tautomer (which builds on [CDK](https://cdk.github.io/)), is a Java library available as a .jar file as well as a [Maven](https://maven.apache.org/) artifact. We need a way to use the Java library in Python and I've found [JPype](https://jpype.readthedocs.io/) to be able to provide what we need.
I haven't thought about fetching the dependencies using Maven itself yet, so the trick right now is to download .jar files. I'm accessing the packages inside the .jar files using JPype.
## Install JPype
Using conda:
```bash
conda install -c conda-forge jpype1
```
If this doesn't work for you, you can find other techniques at [the installation page](https://jpype.readthedocs.io/en/latest/install.html)
## Download the .jar files
Download the .jar files of [Ambit-tautomer](https://sourceforge.net/projects/ambit/files/Ambit2/AMBIT%20applications/tautomers/ambit-tautomers-2.0.0-SNAPSHOT.jar/download) and [CDK](https://github.com/cdk/cdk/releases/tag/cdk-2.3) and place them in the libs/ folder

## Launching the JVM using JPype
A more detailed user guide is available at the [official website](https://jpype.readthedocs.io/en/latest/userguide.html). For our purposes, just do the following: 

Add these imports at the top
```python
# Import the JPype module
import jpype
# Allows Java modules to be imported
import jpype.imports
# Import all standard Java types into the global scope
from jpype.types import *
```

Initialize the JVM and add the .jar files to the classpath

```python
jpype.startJVM(classpath=['libs/ambit-tautomers-2.0.0-SNAPSHOT.jar', 'libs/cdk-2.3.jar'], 
                convertStrings=False)
```

## Accessing Java packages
The native JDK packages should be accessible directly as
```python
from java.lang import System
```
But if you want to access the contents of a .jar file (e.g. the Ambit and CDK .jar files in our case), you need to access it using a method in JPype:

For instance, to import ```ambit2.core.io.FileInputState```, you will have to do
```python
file_input_state = jpype.JPackage('ambit2').core.io.FileInputState
```
## Generating Tautomers
Based on what I understood from their [example code](https://github.com/ideaconsult/apps-ambit/blob/master/tautomers-example/src/main/java/net/idea/example/ambit/tautomers/TautomerWizard.java), the following should be done:

* Create an instance of ```ambit2.tautomers.TautomerManager```, say **tautomerManager**
* Set flags for customization (optional; see next section)
* Create a ```org.openscience.cdk.interfaces.IAtomContainer``` object that will contain the molecule structure, say 'mol'.
    * This could be done using a SMILES string using a ```org.openscience.cdk.smiles.SmilesParser``` instance (see [SmilesParser](https://cdk.github.io/cdk/1.5/docs/api/org/openscience/cdk/smiles/SmilesParser.html)). One could call **SmilesParser#parseSmiles** method from this instance.
* tautomerManager.setStructure(mol)
* Use one of the three algorithms to generate tautomers:
    * Pure combinatorial method using tautomerManager.generateTautomers()
    * Improved combinatorial method using tautomerManager.generateTautomers_ImprovedCombApproach()
    * IA-DFS using tautomerManager.generateTautomersIncrementally()
* These return a ```List<IAtomContainer>``` (a Java object)
### Flags:
TODO
## Problems:
### Creating instances of Java Interfaces
For some purpose, I had to do something along the lines of
```python
silentChemObjectBuilder = jpype.JPackage("org").openscience.cdk.silent.SilentChemObjectBuilder.getInstance()
```
And it threw a ```TypeError: Cannot create Interface instances``` because **SilentChemObjectBuilder** is a Java [Interface](https://docs.oracle.com/javase/tutorial/java/concepts/interface.html) and you cannot instantiate it this way. A [workaround](https://jpype.readthedocs.io/en/latest/userguide.html#case-3-interactive-java) has been suggested in the user guide. Nevertheless, I'm still having some trouble myself.
### Accessing subclasses
```python
GAT = jpype.JPackage("ambit2").tautomers.TautomerConst.GAT
```
throws ```AttributeError: type object 'ambit2.tautomers.TautomerConst' has no attribute 'GAT'```