# BAR/PBSA

Run BAR/PBSA on alchemical decharging trajectories. 
The source explicit solvent simulations should be saved with the listed folder structure where all trajectories and prmtops are stored in folders named by the lambda window. 
A set of example inputs are saved in 1C5X\_inputs.

```
.
└── 1C5X
    └── t1
        ├── complex
        │   ├── 0.000
        │   │   ├── ti.parm7
        │   │   └── ti001.nc
        │   ├── 0.200
        │   │   ├── ti.parm7
        │   │   └── ti001.nc
        │   ├── 0.400
        │   │   ├── ti.parm7
        │   │   └── ti001.nc
        │   ├── 0.600
        │   │   ├── ti.parm7
        │   │   └── ti001.nc
        │   ├── 0.800
        │   │   ├── ti.parm7
        │   │   └── ti001.nc
        │   └── 1.000
        │       ├── ti.parm7
        │       └── ti001.nc
        └── ligands
            ├── 0.000
            │   ├── ti.parm7
            │   └── ti001.nc
            ├── 0.200
            │   ├── ti.parm7
            │   └── ti001.nc
            ├── 0.400
            │   ├── ti.parm7
            │   └── ti001.nc
            ├── 0.600
            │   ├── ti.parm7
            │   └── ti001.nc
            ├── 0.800
            │   ├── ti.parm7
            │   └── ti001.nc
            └── 1.000
                ├── ti.parm7
                └── ti001.nc
```

## Instructions

The BAR/PBSA script automates the 4 stages:

1. Stripping waters and ions from the explicit solvent trajectories. 
Replicate trajectores are concatenated together and RMSD aligned to the first frame.
The ligand is maintained by AMBER mask selection, dummy counter-ions are identified by their lack of charge and maintained.

```
python bar_pbsa.py strip strip_input.yaml
```

2. Writing the PBSA input files with target Radiscale, Protscale, and epsin parameters.

```
python bar_pbsa.py prep prep_input.yaml
```

3. Running the PBSA calculations for neighboring lambdas in parallel. 
The complex and ligand paths are run individually.

```
python bar_pbsa.py run complex_run_input.yaml
python bar_pbsa.py run ligand_run_input.yaml
```

4. Calculating the complete decharging energy with BAR.

```
python bar_pbsa.py calc complex_run_input.yaml
python bar_pbsa.py calc ligand_run_input.yaml
```

