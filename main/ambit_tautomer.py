import jpype
import jpype.imports
from jpype.types import *


jpype.startJVM(classpath=['libs/ambit-tautomers-2.0.0-SNAPSHOT.jar', 'libs/cdk-2.3.jar', 'libs/custom_utils.jar'],
                convertStrings=False)

java = jpype.JPackage("java")
ambit2 = jpype.JPackage("ambit2")
cdk = jpype.JPackage("org").openscience.cdk

tautomerManager = ambit2.tautomers.TautomerManager
#AmbitUtils = jpype.JPackage("kt").AmbitUtilsKt
#AmbitUtils.generateTautomers("CC=O","IA-DFS")
chemObjectBuilder = cdk.silent.SilentChemObjectBuilder.getInstance()
try:
    smilesParser = cdk.smiles.SmilesParser(chemObjectBuilder)
    mol = smilesParser.parseSmiles('CC=O')
except cdk.exception.InvalidSmilesException as e:
    e.printStackTrace()

# It keeps inferring <java org.openscience.cdk.silent.Molecule> as the type of the object
# created above. It needs to be <java org.openscience.cdk.interfaces.IAtomContainer>
# This problem doesn't arise when I try to do this in Kotlin/Java
print(type(mol))

# The following will now throw an exception due to the above bad type inference.
try :
    tautomerManager.setStructure(mol)
    tautomers = tautomerManager.generateTautomersIncrementally()
except java.lang.Exception as e:
    e.printStackTrace()


# Perhaps we could use case labels to decide which algorithm to use.
# (Not) using these constants will not affect code execution. These are used just for reference
#GAT = ambit2.tautomers.TautomerConst.GAT
#print(type(jpype.JClass(ambit2.tautomers.TautomerConst).GAT))
# available algorithms: GAT.Comb_Pure, GAT.Comb_Improved, GAT.Incremental
#algorithm = GAT.Incremental

#rank = TautomerConst.TAUTOMER_RANK
#System.out.println(rank)