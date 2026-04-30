# raspalib — Python API Reference

`raspalib` is a Python binding for the RASPA3 Monte Carlo simulation engine, providing classes for constructing and running grand-canonical and other molecular simulations.

---

## Quick Start

```bash
# Create environment
conda env create -f environment.yml

# Run a single job
python main.py [job.yaml]

# Run batch jobs in parallel
bash run.sh [path/to/jobs] [# of CPUs]
```

---

## Class Overview

| Class | Description |
|---|---|
| [`Atom`](#atom) | Represents a single atom with position, charge, and identity |
| [`Component`](#component) | A molecular component (adsorbate or cation) |
| [`ConnectivityTable`](#connectivitytable) | Molecular connectivity/bonding information |
| [`ForceField`](#forcefield) | Force field parameters and mixing rules |
| [`Framework`](#framework) | Porous framework structure loaded from file |
| [`InputReader`](#inputreader) | Reads simulation configuration from JSON |
| [`IntraMolecularPotentials`](#intramolecularpotentials) | Intra-molecular interaction terms |
| [`Loadings`](#loadings) | Molecule loading counts per component |
| [`MCMoveProbabilities`](#MCMoveProbabilities) | Probabilities for each Monte Carlo move type |
| [`MonteCarlo`](#montecarlo) | Main simulation engine |
| [`PropertyLambdaProbabilityHistogram`](#propertylambdaprobabilityhistogram) | λ-histogram for free energy calculations |
| [`PropertyLoading`](#propertyloading) | Averaged loading statistics over blocks |
| [`PseudoAtom`](#pseudoatom) | Force field atom type definition |
| [`RandomNumber`](#randomnumber) | Seeded random number generator |
| [`RunningEnergy`](#runningenergy) | Accumulates energy contributions |
| [`SimulationBox`](#simulationbox) | Periodic simulation cell (rectangular or triclinic) |
| [`System`](#system) | A complete simulation system |
| [`VDWParameters`](#vdwparameters) | Van der Waals ε/σ parameters |
| [`double3`](#double3) | 3-component float vector |
| [`int3`](#int3) | 3-component integer vector |

---

## API Reference

### Atom

Represents a single particle in the simulation.

**Constructor**

```python
Atom()
Atom(
    position: double3,
    charge: float = 0.0,
    lambda: float = 0.0,
    moleculeId: int = 0,
    type: int = 0,
    componentId: int = 0,
    groupId: bool = False,
    isFractional: bool = False
)
```

**Properties**

| Name | Description |
|---|---|
| `position` | Cartesian coordinates (`double3`) |

---

### Component

A molecular species to be simulated (adsorbate or cation).

**Constructor — from explicit atoms**

```python
Component(
    componentId: int,
    forceField: ForceField,
    componentName: str,
    criticalTemperature: float,
    criticalPressure: float,
    acentricFactor: float,
    definedAtoms: List[Atom],
    connectivitytable: ConnectivityTable = ...,
    intraMolecularPotentials: IntraMolecularPotentials = ...,
    numberOfBlocks: int = 5,
    numberOfLambdaBins: int = 41,
    particleProbabilities: MCMoveProbabilities,
    fugacityCoefficient: Optional[float] = None,
    thermodynamicIntegration: bool = False
)
```

**Constructor — from molecule file**

```python
Component(
    type: Component.Type = Component.Adsorbate,
    componentId: int,
    forceField: ForceField,
    componentName: str,
    fileName: str,
    numberOfBlocks: int = 5,
    numberOfLambdaBins: int = 41,
    particleProbabilities: MCMoveProbabilities = ...,
    fugacityCoefficient: Optional[float] = None,
    thermodynamicIntegration: bool = False
)
```

**Properties**

| Name | Access | Description |
|---|---|---|
| `name` | read-only | Component name string |
| `lambdaGC` | read-only | λ parameter for CFCMC |
| `mc_moves_statistics` | read-only | MC move acceptance statistics |

**Enum: `Component.Type`**

| Value | Description |
|---|---|
| `Component.Adsorbate` | Guest molecule |
| `Component.Cation` | Cation species |

---

### ConnectivityTable

Stores molecular bonding/connectivity information.

```python
ConnectivityTable()
```

---

### ForceField

Defines inter-atomic interaction parameters.

**Constructor — explicit**

```python
ForceField(
    pseudoAtoms: List[PseudoAtom],
    parameters: List[VDWParameters],
    mixingRule: ForceField.MixingRule,
    cutOffFrameworkVDW: float = 12.0,
    cutOffMoleculeVDW: float = 12.0,
    cutOffCoulomb: float = 12.0,
    shifted: bool = True,
    tailCorrections: bool = False,
    useCharge: bool = True
)
```

**Constructor — from file**

```python
ForceField(fileName: str = 'force_field.json')
```

**Properties**

| Name | Access | Description |
|---|---|---|
| `pseudoAtoms` | read-only | List of pseudo-atom types |
| `vdwParameters` | read-only | VDW parameter table |
| `useCharge` | read/write | Enable electrostatics |

**Enum: `ForceField.MixingRule`**

| Value | Description |
|---|---|
| `ForceField.Lorentz_Berthelot` | Lorentz-Berthelot combining rules |

---

### Framework

A porous solid framework structure.

**Constructor — explicit**

```python
Framework(
    frameworkId: int,
    forceField: ForceField,
    componentName: str,
    simulationBox: SimulationBox,
    spaceGroupHallNumber: int,
    definedAtoms: List[Atom],
    numberOfUnitCells: int3
)
```

**Constructor — from file**

```python
Framework(
    frameworkId: int,
    forceField: ForceField,
    componentName: str,
    fileName: Optional[str],
    numberOfUnitCells: int3,
    useChargesFrom: Framework.UseChargesFrom = Framework.PseudoAtoms
)
```

**Properties**

| Name | Access | Description |
|---|---|---|
| `name` | read-only | Framework name |

**Enum: `Framework.UseChargesFrom`**

| Value | Description |
|---|---|
| `Framework.PseudoAtoms` | Use charges from pseudo-atom definitions |

---

### InputReader

Parses a `simulation.json` configuration file.

```python
InputReader(fileName: str = 'simulation.json')
```

**Properties**

| Name | Description |
|---|---|
| `forceField` | Loaded `ForceField` |
| `systems` | List of configured `System` objects |
| `numberOfCycles` | Production cycle count |
| `numberOfInitializationCycles` | Initialization cycle count |
| `numberOfEquilibrationCycles` | Equilibration cycle count |
| `numberOfBlocks` | Statistical block count |
| `printEvery` | Output interval (cycles) |
| `writeEvery` | File write interval (cycles) |
| `writeBinaryRestartEvery` | Restart write interval |
| `optimizeMCMovesEvery` | MC move optimization interval |
| `rescaleWangLandauEvery` | Wang-Landau rescaling interval |

**Enum: `InputReader.SimulationType`**

| Value | Int |
|---|---|
| `MonteCarlo` | 0 |
| `MolecularDynamics` | 2 |
| `Minimization` | 3 |
| `Test` | 4 |
| `Breakthrough` | 5 |
| `MixturePrediction` | 6 |
| `Fitting` | 7 |
| `ParallelTempering` | 8 |

---

### IntraMolecularPotentials

Placeholder for intra-molecular bonded interactions.

```python
IntraMolecularPotentials()
```

---

### Loadings

Tracks the number of molecules per component.

```python
Loadings(numberOfComponents: int)
```

**Methods**

```python
printStatus(component, pressure, unitCells) -> str
printStatus(component, loadings, loadings2, pressure, unitCells) -> str
```

**Properties**

| Name | Description |
|---|---|
| `numberOfMolecules` | Current molecule count |

---

### MCMoveProbabilities

Specifies relative probabilities for each MC move type.

```python
MCMoveProbabilities(
    translationProbability: float = 0.0,
    randomTranslationProbability: float = 0.0,
    rotationProbability: float = 0.0,
    randomRotationProbability: float = 0.0,
    volumeChangeProbability: float = 0.0,
    reinsertionCBMCProbability: float = 0.0,
    swapProbability: float = 0.0,
    swapCBMCProbability: float = 0.0,
    swapCFCMCProbability: float = 0.0,
    widomProbability: float = 0.0,
    gibbsSwapCBMCProbability: float = 0.0,
    # ... additional move types
)
```

**Methods**

```python
join(other: MCMoveProbabilities) -> None   # merge two probability sets
```

---

### MonteCarlo

The primary simulation driver.

**Constructor — explicit**

```python
MonteCarlo(
    numberOfCycles: int,
    numberOfInitializationCycles: int,
    numberOfEquilibrationCycles: int = 0,
    printEvery: int = 5000,
    writeBinaryRestartEvery: int = 5000,
    rescaleWangLandauEvery: int = 1000,
    optimizeMCMovesEvery: int = 100,
    systems: List[System],
    randomSeed: RandomNumber = ...,
    numberOfBlocks: int = 5,
    outputToFiles: bool = False
)
```

**Constructor — from InputReader**

```python
MonteCarlo(inputReader: InputReader)
```

**Methods**

| Method | Description |
|---|---|
| `run()` | Execute full simulation (init + equilibrate + production) |
| `initialize()` | Run initialization phase only |
| `equilibrate()` | Run equilibration phase only |
| `production()` | Run production phase only |
| `cycle()` | Execute a single MC cycle |

**Properties**

| Name | Access | Description |
|---|---|---|
| `systems` | read-only | List of `System` objects |
| `simulationStage` | read/write | Current `SimulationStage` |

**Enum: `MonteCarlo.SimulationStage`**

| Value | Int |
|---|---|
| `Uninitialized` | 0 |
| `Initialization` | 1 |
| `Equilibration` | 2 |
| `Production` | 3 |

---

### PropertyLambdaProbabilityHistogram

Tracks λ-probability histograms for free energy calculations.

```python
PropertyLambdaProbabilityHistogram()
```

**Methods**

```python
normalizedAverageProbabilityHistogram() -> Tuple[List[float], List[float]]
```

**Properties**

| Name | Description |
|---|---|
| `histogram` | Raw histogram data |
| `biasFactor` | Wang-Landau bias factors |

---

### PropertyLoading

Accumulates block-averaged loading statistics.

```python
PropertyLoading(numberOfBlocks: int, numberOfComponents: int)
```

**Methods**

```python
averageLoading() -> Tuple[Loadings, Loadings]
averageLoadingNumberOfMolecules(componentIndex: int) -> Tuple[float, float]
writeAveragesStatistics(components, pressure, unitCells) -> str
```

---

### PseudoAtom

Defines a force field atom type.

```python
PseudoAtom(
    name: str,
    frameworkType: bool,
    mass: float,
    charge: float,
    polarizability: float = 0.0,
    atomicNumber: int,
    printToPDB: bool = True,
    source: str = ''
)
```

---

### RandomNumber

Seeded pseudo-random number generator.

```python
RandomNumber(seed: int = 12)
```

---

### RunningEnergy

Accumulates decomposed energy contributions during simulation.

```python
RunningEnergy()
```

**Properties**

| Name | Description |
|---|---|
| `frameworkMoleculeVDW` | Framework–molecule VDW energy |
| `moleculeMoleculeVDW` | Molecule–molecule VDW energy |

---

### SimulationBox

Defines the periodic simulation cell.

**Constructors**

```python
SimulationBox(a, b, c)                              # orthogonal
SimulationBox(a, b, c, alpha, beta, gamma)          # triclinic angles
SimulationBox(a, b, c, alpha, beta, gamma, type)    # explicit type
SimulationBox(matrix: double3x3, type)              # from cell matrix
```

**Properties**

| Name | Access | Description |
|---|---|---|
| `lengthA` | read-only | Cell length a (Å) |
| `lengthB` | read-only | Cell length b (Å) |
| `lengthC` | read-only | Cell length c (Å) |
| `type` | read/write | `Rectangular` or `Triclinic` |

**Enum: `SimulationBox.Type`**

| Value | Description |
|---|---|
| `SimulationBox.Rectangular` | Orthogonal cell |
| `SimulationBox.Triclinic` | General triclinic cell |

---

### System

A complete simulation system combining framework, components, and conditions.

```python
System(
    systemId: int,
    forceField: ForceField,
    simulationBox: Optional[SimulationBox] = None,
    externalTemperature: float,
    externalPressure: Optional[float] = None,
    heliumVoidFraction: float = 0.0,
    frameworkComponents: Optional[List[Framework]] = [],
    components: List[Component],
    initialPositions: List[List[double3]] = [],
    initialNumberOfMolecules: List[int] = [],
    numberOfBlocks: int = 5,
    systemProbabilities: MCMoveProbabilities = ...,
    sampleMoviesEvery: Optional[int] = None
)
```

**Methods**

| Method | Description |
|---|---|
| `computeTotalEnergies()` | Returns current `RunningEnergy` |
| `frameworkMass()` | Returns framework mass (g/mol) or `None` |
| `writeMCMoveStatistics()` | Returns formatted move statistics string |

**Properties**

| Name | Access | Description |
|---|---|---|
| `components` | read-only | List of `Component` objects |
| `loadings` | read-only | Current `Loadings` |
| `averageLoadings` | read-only | Block-averaged `PropertyLoading` |
| `inputPressure` | read-only | Specified external pressure |
| `atomData` | read/write | Raw atom position/type array |

---

### VDWParameters

Lennard-Jones ε and σ parameters for a pair interaction.

```python
VDWParameters(epsilon: float, sigma: float)
```

---

### double3

A 3-component vector of `float` values.

```python
double3(x: float = 0.0, y: float = 0.0, z: float = 0.0)
```

**Properties:** `x`, `y`, `z`

---

### int3

A 3-component vector of `int` values.

```python
int3(x: int = 0, y: int = 0, z: int = 0)
```

**Properties:** `x`, `y`, `z`

---

## Module Location

```
/miniconda3/envs/raspa3/lib/python3.9/site-packages/raspalib.cpython-39-x86_64-linux-gnu.so
```
