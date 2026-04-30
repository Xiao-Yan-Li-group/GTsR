conda env create -f environment.yml
python main.py [job.yaml]
bash run.sh [path/to/jobs] [# of CPUs]

Help on module raspalib:

NAME
    raspalib

CLASSES
    pybind11_builtins.pybind11_object(builtins.object)
        Atom
        Component
        ConnectivityTable
        ForceField
        Framework
        InputReader
        IntraMolecularPotentials
        Loadings
        MCMoveProbabilities
        MonteCarlo
        PropertyLambdaProbabilityHistogram
        PropertyLoading
        PseudoAtom
        RandomNumber
        RunningEnergy
        SimulationBox
        System
        VDWParameters
        double3
        int3

    class Atom(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      Atom
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(*args, **kwargs)
     |      Overloaded function.
     |
     |      1. __init__(self: raspalib.Atom) -> None
     |
     |      2. __init__(self: raspalib.Atom, position: raspalib.double3, charge: float = 0.0, lambda: float = 0.0, moleculeId: int = 0, type: int = 0, componentId: int = 0, groupId: bool = 0, isFractional: bool = 0) -> None
     |
     |  __repr__(...)
     |      __repr__(self: raspalib.Atom) -> str
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  position
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class Component(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      Component
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(*args, **kwargs)
     |      Overloaded function.
     |
     |      1. __init__(self: raspalib.Component, componentId: int, forceField: raspalib.ForceField, componentName: str, criticalTemperature: float, criticalPressure: float, acentricFactor: float, definedAtoms: List[raspalib.Atom], connectivitytable: raspalib.ConnectivityTable = <raspalib.ConnectivityTable object at 0x79e4b9e41d30>, intraMolecularPotentials: raspalib.IntraMolecularPotentials = <raspalib.IntraMolecularPotentials object at 0x79e4b9e41d70>, numberOfBlocks: int = 5, numberOfLambdaBins: int = 41, particleProbabilities: raspalib.MCMoveProbabilities, fugacityCoefficient: Optional[float] = None, thermodynamicIntegration: bool = False) -> None
     |
     |      2. __init__(self: raspalib.Component, type: raspalib.Component.Type = <Type.Adsorbate: 0>, componentId: int, forceField: raspalib.ForceField, componentName: str, fileName: str, numberOfBlocks: int = 5, numberOfLambdaBins: int = 41, particleProbabilities: raspalib.MCMoveProbabilities = <raspalib.MCMoveProbabilities object at 0x79e4b9e41db0>, fugacityCoefficient: Optional[float] = None, thermodynamicIntegration: bool = False) -> None
     |
     |  __repr__(...)
     |      __repr__(self: raspalib.Component) -> str
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  lambdaGC
     |
     |  mc_moves_statistics
     |
     |  name
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  Adsorbate = <Type.Adsorbate: 0>
     |
     |  Cation = <Type.Cation: 1>
     |
     |  Type = <class 'raspalib.Component.Type'>
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class ConnectivityTable(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      ConnectivityTable
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.ConnectivityTable) -> None
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class ForceField(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      ForceField
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(*args, **kwargs)
     |      Overloaded function.
     |
     |      1. __init__(self: raspalib.ForceField, pseudoAtoms: List[raspalib.PseudoAtom], parameters: List[raspalib.VDWParameters], mixingRule: ForceField@forcefield::MixingRule, cutOffFrameworkVDW: float = 12.0, cutOffMoleculeVDW: float = 12.0, cutOffCoulomb: float = 12.0, shifted: bool = True, tailCorrections: bool = False, useCharge: bool = True) -> None
     |
     |      2. __init__(self: raspalib.ForceField, fileName: str = 'force_field.json') -> None
     |
     |  __repr__(...)
     |      __repr__(self: raspalib.ForceField) -> str
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  pseudoAtoms
     |
     |  vdwParameters
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  useCharge
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  Lorentz_Berthelot = <MixingRule.Lorentz_Berthelot: 0>
     |
     |  MixingRule = <class 'raspalib.ForceField.MixingRule'>
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class Framework(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      Framework
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(*args, **kwargs)
     |      Overloaded function.
     |
     |      1. __init__(self: raspalib.Framework, frameworkId: int, forceField: raspalib.ForceField, componentName: str, simulationBox: raspalib.SimulationBox, spaceGroupHallNumber: int, definedAtoms: List[raspalib.Atom], numberOfUnitCells: raspalib.int3) -> None
     |
     |      2. __init__(self: raspalib.Framework, frameworkId: int, forceField: raspalib.ForceField, componentName: str, fileName: Optional[str], numberOfUnitCells: raspalib.int3, useChargesFrom: raspalib.Framework.UseChargesFrom = <UseChargesFrom.PseudoAtoms: 0>) -> None
     |
     |  __repr__(...)
     |      __repr__(self: raspalib.Framework) -> str
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  name
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  UseChargesFrom = <class 'raspalib.Framework.UseChargesFrom'>
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class InputReader(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      InputReader
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.InputReader, fileName: str = 'simulation.json') -> None
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  forceField
     |
     |  numberOfBlocks
     |
     |  numberOfCycles
     |
     |  numberOfEquilibrationCycles
     |
     |  numberOfInitializationCycles
     |
     |  optimizeMCMovesEvery
     |
     |  printEvery
     |
     |  rescaleWangLandauEvery
     |
     |  systems
     |
     |  writeBinaryRestartEvery
     |
     |  writeEvery
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  Breakthrough = <SimulationType.Breakthrough: 5>
     |
     |  Fitting = <SimulationType.Fitting: 7>
     |
     |  Minimization = <SimulationType.Minimization: 3>
     |
     |  MixturePrediction = <SimulationType.MixturePrediction: 6>
     |
     |  MolecularDynamics = <SimulationType.MolecularDynamics: 2>
     |
     |  MonteCarlo = <SimulationType.MonteCarlo: 0>
     |
     |  MonteCarloTransitionMatrix = <SimulationType.MonteCarloTransitionMatri...
     |
     |  ParallelTempering = <SimulationType.ParallelTempering: 8>
     |
     |  SimulationType = <class 'raspalib.InputReader.SimulationType'>
     |
     |  Test = <SimulationType.Test: 4>
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class IntraMolecularPotentials(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      IntraMolecularPotentials
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.IntraMolecularPotentials) -> None
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class Loadings(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      Loadings
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.Loadings, arg0: int) -> None
     |
     |  printStatus(...)
     |      printStatus(*args, **kwargs)
     |      Overloaded function.
     |
     |      1. printStatus(self: raspalib.Loadings, arg0: raspalib.Component, arg1: Optional[float], arg2: Optional[raspalib.int3]) -> str
     |
     |      2. printStatus(self: raspalib.Loadings, arg0: raspalib.Component, arg1: raspalib.Loadings, arg2: raspalib.Loadings, arg3: Optional[float], arg4: Optional[raspalib.int3]) -> str
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  numberOfMolecules
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class MCMoveProbabilities(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      MCMoveProbabilities
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.MCMoveProbabilities, translationProbability: float = 0.0, randomTranslationProbability: float = 0.0, rotationProbability: float = 0.0, randomRotationProbability: float = 0.0, volumeChangeProbability: float = 0.0, reinsertionCBMCProbability: float = 0.0, identityChangeProbability: float = 0.0, swapProbability: float = 0.0, swapCBMCProbability: float = 0.0, swapCFCMCProbability: float = 0.0, swapCBCFCMCProbability: float = 0.0, gibbsVolumeChangeProbability: float = 0.0, gibbsSwapCBMCProbability: float = 0.0, gibbsSwapCFCMCProbability: float = 0.0, widomProbability: float = 0.0, widomCFCMCProbability: float = 0.0, widomCBCFCMCProbability: float = 0.0, parallelTemperingProbability: float = 0.0, hybridMCProbability: float = 0.0) -> None
     |
     |  join(...)
     |      join(self: raspalib.MCMoveProbabilities, arg0: raspalib.MCMoveProbabilities) -> None
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class MonteCarlo(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      MonteCarlo
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(*args, **kwargs)
     |      Overloaded function.
     |
     |      1. __init__(self: raspalib.MonteCarlo, numberOfCycles: int, numberOfInitializationCycles: int, numberOfEquilibrationCycles: int = 0, printEvery: int = 5000, writeBinaryRestartEvery: int = 5000, rescaleWangLandauEvery: int = 1000, optimizeMCMovesEvery: int = 100, systems: List[raspalib.System], randomSeed: raspalib.RandomNumber = <raspalib.RandomNumber object at 0x79e4b9e48a30>, numberOfBlocks: int = 5, outputToFiles: bool = False) -> None
     |
     |      2. __init__(self: raspalib.MonteCarlo, inputReader: raspalib.InputReader) -> None
     |
     |  cycle(...)
     |      cycle(self: raspalib.MonteCarlo) -> None
     |
     |  equilibrate(...)
     |      equilibrate(self: raspalib.MonteCarlo) -> None
     |
     |  initialize(...)
     |      initialize(self: raspalib.MonteCarlo) -> None
     |
     |  production(...)
     |      production(self: raspalib.MonteCarlo) -> None
     |
     |  run(...)
     |      run(self: raspalib.MonteCarlo) -> None
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  systems
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  simulationStage
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  Equilibration = <SimulationStage.Equilibration: 2>
     |
     |  Initialization = <SimulationStage.Initialization: 1>
     |
     |  Production = <SimulationStage.Production: 3>
     |
     |  SimulationStage = <class 'raspalib.MonteCarlo.SimulationStage'>
     |
     |  Uninitialized = <SimulationStage.Uninitialized: 0>
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class PropertyLambdaProbabilityHistogram(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      PropertyLambdaProbabilityHistogram
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.PropertyLambdaProbabilityHistogram) -> None
     |
     |  normalizedAverageProbabilityHistogram(...)
     |      normalizedAverageProbabilityHistogram(self: raspalib.PropertyLambdaProbabilityHistogram) -> Tuple[List[float], List[float]]
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  biasFactor
     |
     |  histogram
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class PropertyLoading(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      PropertyLoading
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.PropertyLoading, arg0: int, arg1: int) -> None
     |
     |  __repr__(...)
     |      __repr__(self: raspalib.PropertyLoading) -> str
     |
     |  averageLoading(...)
     |      averageLoading(self: raspalib.PropertyLoading) -> Tuple[raspalib.Loadings, raspalib.Loadings]
     |
     |  averageLoadingNumberOfMolecules(...)
     |      averageLoadingNumberOfMolecules(self: raspalib.PropertyLoading, arg0: int) -> Tuple[float, float]
     |
     |  writeAveragesStatistics(...)
     |      writeAveragesStatistics(self: raspalib.PropertyLoading, arg0: List[raspalib.Component], arg1: Optional[float], arg2: Optional[raspalib.int3]) -> str
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class PseudoAtom(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      PseudoAtom
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.PseudoAtom, name: str, frameworkType: bool, mass: float, charge: float, polarizability: float = 0.0, atomicNumber: int, printToPDB: bool = True, source: str = '') -> None
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class RandomNumber(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      RandomNumber
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.RandomNumber, seed: int = 12) -> None
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class RunningEnergy(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      RunningEnergy
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.RunningEnergy) -> None
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  frameworkMoleculeVDW
     |
     |  moleculeMoleculeVDW
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class SimulationBox(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      SimulationBox
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(*args, **kwargs)
     |      Overloaded function.
     |
     |      1. __init__(self: raspalib.SimulationBox, a: float, b: float, c: float) -> None
     |
     |      2. __init__(self: raspalib.SimulationBox, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> None
     |
     |      3. __init__(self: raspalib.SimulationBox, a: float, b: float, c: float, alpha: float, beta: float, gamma: float, type: raspalib.SimulationBox.Type) -> None
     |
     |      4. __init__(self: raspalib.SimulationBox, arg0: double3x3@double3x3, arg1: raspalib.SimulationBox.Type) -> None
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  lengthA
     |
     |  lengthB
     |
     |  lengthC
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  type
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  Rectangular = <Type.Rectangular: 0>
     |
     |  Triclinic = <Type.Triclinic: 1>
     |
     |  Type = <class 'raspalib.SimulationBox.Type'>
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class System(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      System
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.System, systemId: int, forceField: raspalib.ForceField, simulationBox: Optional[raspalib.SimulationBox] = None, externalTemperature: float, externalPressure: Optional[float] = None, heliumVoidFraction: float = 0.0, frameworkComponents: Optional[raspalib.Framework] = [], components: List[raspalib.Component], initialPositions: List[List[raspalib.double3]] = [], initialNumberOfMolecules: List[int] = [], numberOfBlocks: int = 5, systemProbabilities: raspalib.MCMoveProbabilities = <raspalib.MCMoveProbabilities object at 0x79e4b9e43cb0>, sampleMoviesEvery: Optional[int] = None) -> None
     |
     |  __repr__(...)
     |      __repr__(self: raspalib.System) -> str
     |
     |  computeTotalEnergies(...)
     |      computeTotalEnergies(self: raspalib.System) -> raspalib.RunningEnergy
     |
     |  frameworkMass(...)
     |      frameworkMass(self: raspalib.System) -> Optional[float]
     |
     |  writeMCMoveStatistics(...)
     |      writeMCMoveStatistics(self: raspalib.System) -> str
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  averageLoadings
     |
     |  components
     |
     |  inputPressure
     |
     |  loadings
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  atomData
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class VDWParameters(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      VDWParameters
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.VDWParameters, epsilon: float, sigma: float) -> None
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class double3(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      double3
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.double3, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> None
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  x
     |
     |  y
     |
     |  z
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

    class int3(pybind11_builtins.pybind11_object)
     |  Method resolution order:
     |      int3
     |      pybind11_builtins.pybind11_object
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(...)
     |      __init__(self: raspalib.int3, x: int = 0, y: int = 0, z: int = 0) -> None
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  x
     |
     |  y
     |
     |  z
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from pybind11_builtins.pybind11_object:
     |
     |  __new__(*args, **kwargs) from pybind11_builtins.pybind11_type
     |      Create and return a new object.  See help(type) for accurate signature.

FILE
    /miniconda3/envs/raspa3/lib/python3.9/site-packages/raspalib.cpython-39-x86_64-linux-gnu.so
