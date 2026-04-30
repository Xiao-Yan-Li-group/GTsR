import raspalib


def setting(system,
            n_cycles=20000,
            n_cycles_init=10000,
            n_cycles_eq=0,
            print_every=5000,
            write_every=5000,
            recale_every=1000,
            optimize_every=100,
            random_seed=None,
            n_blocks=5,
            save=True
            ):
    
    if random_seed is None:
        random_seed = raspalib.RandomNumber(12)

    return raspalib.MonteCarlo(
                    numberOfCycles=n_cycles,
                    numberOfInitializationCycles=n_cycles_init,
                    numberOfEquilibrationCycles = n_cycles_eq,
                    printEvery=print_every,
                    writeBinaryRestartEvery=write_every,
                    rescaleWangLandauEvery=recale_every,
                    optimizeMCMovesEvery=optimize_every,
                    systems=[system],
                    randomSeed=random_seed,
                    numberOfBlocks=n_blocks,
                    outputToFiles=save
                )


def job(mc):
    mc.run()