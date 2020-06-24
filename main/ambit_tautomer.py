import os
import jpype
import jpype.imports
from jpype import JImplements, JOverride
from jpype.types import *

jpype.startJVM(classpath=['libs/ambit-tautomers-2.0.0-SNAPSHOT.jar', 'libs/cdk-2.3.jar'], 
                convertStrings=False)
from java.lang import System

# I was printing the classpath for testing
#print(System.getProperty("java.class.path"))

ambit2 = jpype.JPackage("ambit2")
cdk = jpype.JPackage("org").openscience.cdk

tautomerManager = ambit2.tautomers.TautomerManager()

# If you're wondering why I did this, see Using_Ambit.md (Java Interfaces can't be instantiated)
@JImplements(cdk.interfaces.IChemObjectBuilder)
class PythonCDKObjectBuilder:
    def __init__(self):
        # I don't know if anything else needs to be done but this compiles fine
        pass
    # Must override Interface method "newInstance"
    @JOverride
    def newInstance(object):
        return self
objectBuilder = PythonCDKObjectBuilder()
smilesParser = cdk.smiles.SmilesParser(objectBuilder)

# TODO: fix this
# The following line throws the seemingly unrelated error
# TypeError: newInstance() takes 1 positional argument but 3 were given
mol = smilesParser.parseSmiles('CC=O')

#tautomerManager.setStructure(mol)
#tautomerManager.generateTautomersIncrementaly

# Perhaps we could use case labels to decide which algorithm to use.
# (Not) using these constants will not affect code execution. These are used just for reference
#GAT = ambit2.tautomers.TautomerConst.GAT
#print(type(jpype.JClass(ambit2.tautomers.TautomerConst).GAT))
# available algorithms: GAT.Comb_Pure, GAT.Comb_Improved, GAT.Incremental
#algorithm = GAT.Incremental

#rank = TautomerConst.TAUTOMER_RANK
#System.out.println(rank)