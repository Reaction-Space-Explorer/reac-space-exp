# Using Ambit in Python
I'm accessing the contents of the .jar files using JPype
## Install JPype
Using conda:
```bash
conda install -c conda-forge jpype1
```
If this doesn't work for you, you can find other techniques at [the installation page](https://jpype.readthedocs.io/en/latest/install.html)
## Download the .jar files
Download the .jar files of [Ambit-tautomer](https://sourceforge.net/projects/ambit/files/Ambit2/AMBIT%20applications/tautomers/ambit-tautomers-2.0.0-SNAPSHOT.jar/download) and [CDK](https://github.com/cdk/cdk/releases/tag/cdk-2.3) and place them in the libs/ folder

## Launching the JVM using JPype
Add these imports at the top
```python
import jpype
import jpype.imports
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