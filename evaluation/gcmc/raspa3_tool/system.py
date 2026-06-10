import raspalib
from pathlib import Path


def framework(id,
              cif_file,
              uc,
              ff):
    cif_file = Path(cif_file).resolve()
    return raspalib.Framework(
        frameworkId=id,
        forceField=ff,
        componentName=cif_file.stem,
        fileName=str(cif_file),
        numberOfUnitCells=raspalib.int3(uc[0], uc[1], uc[2]),
        useChargesFrom=raspalib.Framework.UseChargesFrom.CIF_File
    )

def adsorbate(id,
              ff,
              file,
              probs,
              n_blocks=5,
              n_bins=41,
              f=1.0,
              thermodynamicIntegration=False):
    gas = raspalib.Component(
        type=raspalib.Component.Adsorbate,
        componentId=id,
        forceField=ff,
        componentName=Path(file).stem,
        fileName=file,
        numberOfBlocks=n_blocks,
        numberOfLambdaBins=n_bins,
        particleProbabilities=probs,
        fugacityCoefficient=f,
        thermodynamicIntegration=thermodynamicIntegration
        )
    return gas

def probs(
    translationProbability=0.0,
    randomTranslationProbability=0.0,
    rotationProbability=0.0,
    randomRotationProbability=0.0,
    volumeChangeProbability=0.0,
    reinsertionCBMCProbability=0.0,
    identityChangeProbability=0.0,
    swapProbability=0.0,
    swapCBMCProbability=0.0,
    swapCFCMCProbability=0.0,
    swapCBCFCMCProbability=0.0,
    gibbsVolumeChangeProbability=0.0,
    gibbsSwapCBMCProbability=0.0,
    gibbsSwapCFCMCProbability=0.0,
    widomProbability=0.0,
    widomCFCMCProbability=0.0,
    widomCBCFCMCProbability=0.0,
    parallelTemperingProbability=0.0,
    hybridMCProbability=0.0,
):
    return raspalib.MCMoveProbabilities(
        translationProbability=translationProbability,
        randomTranslationProbability=randomTranslationProbability,
        rotationProbability=rotationProbability,
        randomRotationProbability=randomRotationProbability,
        volumeChangeProbability=volumeChangeProbability,
        reinsertionCBMCProbability=reinsertionCBMCProbability,
        identityChangeProbability=identityChangeProbability,
        swapProbability=swapProbability,
        swapCBMCProbability=swapCBMCProbability,
        swapCFCMCProbability=swapCFCMCProbability,
        swapCBCFCMCProbability=swapCBCFCMCProbability,
        gibbsVolumeChangeProbability=gibbsVolumeChangeProbability,
        gibbsSwapCBMCProbability=gibbsSwapCBMCProbability,
        gibbsSwapCFCMCProbability=gibbsSwapCFCMCProbability,
        widomProbability=widomProbability,
        widomCFCMCProbability=widomCFCMCProbability,
        widomCBCFCMCProbability=widomCBCFCMCProbability,
        parallelTemperingProbability=parallelTemperingProbability,
        hybridMCProbability=hybridMCProbability,
    )

def setting(    
            id,
            ff,
            temperature,
            pressure,
            framework,
            components,
            n_adsorbates_init,
            sysMC,
            simulationBox=None,
            n_blocks=5,
            pos_init=None,
            VF_He=0.0,
            movies_every=None,
            ):
    
    if pos_init is None:
        pos_init = []
        
    return raspalib.System(
            systemId=id,
            forceField=ff,
            externalTemperature=temperature,
            externalPressure=pressure,
            frameworkComponents=framework,
            components=components,
            initialNumberOfMolecules=n_adsorbates_init,
            simulationBox=simulationBox,
            numberOfBlocks=n_blocks,
            initialPositions=pos_init,
            heliumVoidFraction=VF_He,
            sampleMoviesEvery=movies_every,
            systemProbabilities=sysMC
            )
