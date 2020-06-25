import jpype
import jpype.imports
from jpype.types import *


jpype.startJVM(classpath=['libs/ambit-tautomers-2.0.0-SNAPSHOT.jar', 'libs/cdk-2.3.jar'],
                convertStrings=False)

java = jpype.JPackage("java")
ambit2 = jpype.JPackage("ambit2")
cdk = jpype.JPackage("org").openscience.cdk

def generateTautomers(smiles, mode="IA-DFS"):
    """
    Generate the list of possible tautomers for a given molecule

    Keyword arguments:
    smiles -- the SMILES string of the molecule
    mode -- The generation algorithm to use (default: "IA-DFS")

    Available  algorithms include:
        "CM" for Simple Combinatorial,
        "CMI" for Improved Combinatorial,
        "IA-DFS" for Incremental Algorithm - Depth First Search
        "combined" for a combination of CMI and IA-DFS
    See https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.201200133 for a detailed discussion
        of the algorithms

    Returns: a list of SMILES strings of the possible tautomers
    """
    tautomerManager = ambit2.tautomers.TautomerManager()
    silentChemObjectBuilder = cdk.silent.SilentChemObjectBuilder.getInstance()
    try:
        smilesParser = cdk.smiles.SmilesParser(silentChemObjectBuilder)
        mol = cdk.interfaces.IAtomContainer
        mol = smilesParser.parseSmiles('CC=O')
    except cdk.exception.CDKException as e:
        e.printStackTrace()

    tautomerManager.setStructure(mol)
    if mode == "IA-DFS":
        tautomers = tautomerManager.generateTautomersIncrementaly()
    elif mode == "combined":
        tautomers = tautomerManager.generateTautomersCombinedApproach()
    elif mode == "CMI":
        tautomers = tautomerManager.generateTautomers_ImprovedCombApproach()
    elif mode == "CM":
        tautomers = tautomerManager.generateTautomers()
    smilesGenerator = cdk.smiles.SmilesGenerator(True)
    return [smilesGenerator.createSMILES(taut) for taut in tautomers]