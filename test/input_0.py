database(
    thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0'], # it can't predict reactions without them
 #   seedMechanisms = [],
 #    kineticsDepositories=['training'],
 #   kineticsFamilies='default',
 #   kineticsEstimator='rate rules',
    reactionLibraries = [] # 'Glarborg/C3' had forbidden Species C20 and threw an exception
)

## isomeric smiles for glucose: C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O
## canonical smiles for glucose: 'C(C1C(C(C(C(O1)O)O)O)O)O'
species(
    label='glycine',
    reactive=True,
    structure=SMILES('NCC(=O)O'),
)

species(
    label='HCN',
    reactive=True,
    structure=SMILES('C#N')
)
species(
    label='glucose',
    reactive=True,
    structure=SMILES('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O')
)

simpleReactor(
    temperature=(800,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        'glucose': 0.4,
        'glycine':0.5,
        'HCN':0.1
    },
    terminationTime=(1e0,'s'),
    terminationConversion={
        'glycine':0.99, # terminate the reaction upon 99% conversion of glycine
        'glucose':0.99
    },
    sensitivity=('HCN','glycine','glucose')
)

## Absolute and relative tolerance(s)
## I don't understand what exactly these do
simulator(
    atol=1e-16,
    rtol=1e-8,
)

#This section must be defined
model(
    # fraction of charcteristic flux an edge species must have (otherwise it'll be pruned)
    toleranceKeepInEdge=0.01,
    # fraction of some "characteristic flux" (rate) an edge species must have for it to be moved to the core
    toleranceMoveToCore=0.1,
    # this si something I still don't understand
    toleranceInterruptSimulation=1,
    # parameters for pruning edge species to make sure they don't exhaust the memory space
    maximumEdgeSpecies=100000,
    minSpeciesExistIterationsForPrune=4,
    filterReactions=True,
)

options(
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=False,
)