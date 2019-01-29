# Notes for Simulations

* restrict attention to 2 wave designs


## Simulations without covariates

* Optimal assignment only works without covariates (or with independent priors)
* Optimal assignment becomes infeasible when second wave is large (or strata in second wave are large)
* Structure for current simulation code without covariates:
		1 ReadAllData() loads data (with and without covariates), and prints figures if prompted
		2 DesignTable() applies Simulations (using ExpectedRegret()) to each theta value, 
			and then produces printable table using PrintRegretTable()
		3 ExpectedRegret() picks a set of assignment methods, and then runs (in parallel across replication draws)
			 the SimulateTWaveDesign() function (or SimulateConventionalDesign() for 1 wave version; but that should maybe be subsumed)
		4 SimulateTWaveDesign() assigns each treatment with same frequency in 1st wave, and then uses Dtchoice 
			(with a methods argument) to get treatment assignment
		5 Dtchoice() applies one of a number of methods to pick treatment vector in the given wave.
* For a fair comparison, SimulateTWaveDesign() permutes theta for each iteration
* modified Thompson method just avoids repeats, but does not average

## Simulations with covariates