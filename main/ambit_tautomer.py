import os
import jpype
import jpype.imports
from jpype.types import *
from rdkit.Chem import MolFromSmiles, MolToSmiles

# Note: this will work if your working directory is reac-space-exp/main
# If your working directory is reac-space-exp/ instead, you need to get rid of the ".."
# Prefer doing "libs/whatever.jar" instead of the entire os.path.join thing
jpype.startJVM(classpath=[os.path.join('..','libs/ambit-tautomers-2.0.0-SNAPSHOT.jar'),
        os.path.join('..','libs/cdk-2.3.jar')], convertStrings=True)

from java.lang import System

java = jpype.JPackage("java")
ambit2 = jpype.JPackage("ambit2")
cdk = jpype.JPackage("org").openscience.cdk

tautomerManager = JClass('ambit2.tautomers.TautomerManager')()
silentChemObjectBuilder = JClass('org.openscience.cdk.silent.SilentChemObjectBuilder').getInstance()


def smilesToMolecule(smiles):
    """
    Convert a SMILES string to a CDK Molecule object.

    Returns: the Molecule object
    """
    mol = None
    try:
        smilesParser = cdk.smiles.SmilesParser(silentChemObjectBuilder)
        mol = smilesParser.parseSmiles(smiles)
    except cdk.exception.InvalidSmilesException as e:
        System.err.println('An error occured while parsing the SMILES')
        e.printStackTrace()
    return mol


def smiles_from_molecule(molecule):
    """
    Parse a CDK object into a SMILES String

    Returns: an *RDKit* canonical SMILES
    """
    smi_flavor = cdk.smiles.SmiFlavor
    smilesGenerator = cdk.smiles.SmilesGenerator(True)
    smiles = smilesGenerator.createSMILES(molecule)
    # Note: this method may not create canonical SMILES prefer using RDKit for canonicalization
    # CDK does have a method to create canonical SMILES but I can't call it using JPype (or even using Java/Kotlin themselves)
    mol = MolFromSmiles(smiles)
    smiles = MolToSmiles(mol)
    return smiles


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
        "combined" doesn't work is intended for a combination of CMI and IA-DFS
    See https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.201200133 for a detailed discussion
        of the algorithms

    Returns: a list of SMILES strings of the possible tautomers
    """
    mol = smilesToMolecule(smiles)
    tautomerManager.setStructure(mol)
    if mode == "IA-DFS":
        tautomers = tautomerManager.generateTautomersIncrementaly()
    # Note: "combined" doesn't produce anything
    elif mode == "combined":
        print("WARNING: Combined approach is not fully complete yet, it may produce nothing")
        tautomers = tautomerManager.generateTautomersCombinedApproach()
    elif mode == "CMI":
        tautomers = tautomerManager.generateTautomers_ImprovedCombApproach()
    elif mode == "CM":
        tautomers = tautomerManager.generateTautomers()
    else:
        raise NameError("Invalid generation algorithm mode: {0}, please correct typos".format(mode))
    smiles_to_return = [smiles_from_molecule(taut) for taut in tautomers]
    # In case the post filter removes the original molecule
    # TODO: Carbons with two double bonds (and other substrucutres) should be forbidden
    # at the time of reaction itself then this shouldn't happen
    if smiles_from_molecule(mol) not in smiles_to_return:
        smiles_to_return.append(smiles_from_molecule(mol))
    return tuple(smiles_to_return)


def setNumBackTracks(num):
    """
    Sets the maximum number of back tracks the IA-DFS algorithm is allowed to perform.
    """
    tautomerManager.maxNumOfBackTracks = num


def setMaxTautomerRegistrations(num):
    tautomerManager.maxNumOfTautomerRegistrations = num


def maxSubCombinations(num):
    """
    Maximum number of sub-combinations the improved combinatorial methods is allowed to work with
    See their paper for details on what a sub-combination is.
    """
    tautomerManager.maxNumOfSubCombiations = num


def toCalculateCACTVSEnergyRank(flag):
    tautomerManager.FlagCalculateCACTVSEnergyRank = flag


def use13Rules(flag):
    """
    Whether or not to use 1,3 shift tautomer rules

    Keyword arguments:
    flag: True or False

    """
    tautomerManager.getKnowledgeBase().FlagUse13Shifts = flag


def use15Rules(flag):
    tautomerManager.getKnowledgeBase().FlagUse15Shifts = flag


def use17Rules(flag):
    tautomerManager.getKnowledgeBase().FlagUse17Shifts = flag


#def selectRules(args):
    # TODO:
    #tautomerManager.getRuleSelector().setSelectionMode(RSM.valueOf(args))

def setRuleNumberLimit(limit):
    tautomerManager.getRuleSelector().setRuleNumberLimit(limit)


def useDuplicationIsomorphismCheck(flag):
    # TODO: see what this even does
    tautomerManager.tautomerFilter.setFlagApplyDuplicationCheckIsomorphism(flag)


def useDuplicationCheckInChI(flag):
    # TODO: same as above.
    tautomerManager.tautomerFilter.setFlagApplyDuplicationCheckInChI(flag)

#jpype.shutdownJVM()