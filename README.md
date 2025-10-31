# pydusty
A generic python wrapper for performing MCMC analysis with the SED modeling code `dusty`.

This code is based on Scott Adam's rewrite of Chris Kochanek's MCMC wrapper around `dusty`, and my subsequent rewrite of Scott's wrapper using `emcee`. The original codes are available under `pydusty/legacy`.

## Installation
1. Clone the repo using `git clone`
2. `cd pydusty`
3. `pip install -e .`

The dustyv2 source code files are under `data/dusty_files`. You may need to recreate the executable file with `gfortran dusty.f -o dusty`. Please see the dusty github repository (esp. https://github.com/ivezic/dusty/issues/7 issue) for more details.

To test your install, copy `data/lums_WISE_J175749.76-075314.9_ebv1.18_nodisterr.dat` to your favorite directory on your computer (hereafter labelled as `<pydusty_data_path>`).<br>
Make a directory `<pydusty_data_path>/dusty_run`. <br>
From the pydusty directory, run `python examples/run_dusty_mcmc.py <pydusty_data_path>/lums_WISE_J175749.76-075314.9_ebv1.18_nodisterr.dat <pydusty_data_path>/dusty_run --loglevel INFO --nprocesses 4`. <br>
The full code takes a couple hours to run, but if it runs without interruptions for ~5 minutes, your install has worked correctly.

If you run into array limit issues with your dusty executable you can update the array bounds like so (I update the limits to 250 but you can do more...):
```
cd pydusty/data/dusty_files
cp dusty.f dusty_original.f
sed -i '' 's/(3,200)/(3,250)/g' dusty.f
```
and recompile like so:
```
gfortran -std=legacy -fallow-argument-mismatch -fno-automatic -fno-align-commons -O2 -o dusty dusty.f
```

(Better tests will follow)

## Usage
An example usage file is given under `examples/run_dusty_mcmc.py (--h)`. Detailed documentation of the code will follow soon.

The following functionalities present in Scott's `dusty/' have not yet been implemented in this version
1. Using observed expansion velocities to constrain dust shell radii
2. Using observed dL/dt as an additional constraint.

It should be straightforward to implement these, I just haven;t encountered data where these are necessary. <br>
Happy to add these in for compelling science cases, or even better, accept contributions here.
