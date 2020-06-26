import jpype
import jpype.imports
from jpype.types import *

# Note: this will work if your working directory is reac-space-exp
# If you're in reac-space-exp/main instead, you need to add a ".." before the paths
# Prefer doing os.path.join("..", "libs/whatever.jar")
jpype.startJVM(classpath=['libs/ambit-tautomers-2.0.0-SNAPSHOT.jar', 'libs/cdk-2.3.jar'],
                convertStrings=False)

java = jpype.JPackage("java")
ambit2 = jpype.JPackage("ambit2")
cdk = jpype.JPackage("org").openscience.cdk

tautomerManager = ambit2.tautomers.TautomerManager()


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
    silentChemObjectBuilder = cdk.silent.SilentChemObjectBuilder.getInstance()
    try:
        smilesParser = cdk.smiles.SmilesParser(silentChemObjectBuilder)
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


def setNumBackTracks(num):
    tautomerManager.maxNumOfBackTracks = num


def setMaxTautomerRegistrations(num):
    tautomerManager.maxNumOfTautomerRegistrations = num


def maxSubCombinations(num):
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