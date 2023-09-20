# CG_IDP (Coarse-grained model of IDPs)

The code uses a hydropathy parameter called M3, which was developed by Tesai et al. in 2021 ([DOI Link](https://doi.org/10.1073/pnas.2111696118)). 
To run the code for 33 IDPS (listed in idp_bank.csv), please follow the following procedure:
1. Run master_IDPs.py, and when prompted, specify the index of the IDP sequence to run. This will create a directory with the name of the IDP.
2. cd to the directory and run master.py. This will generate folders with epsilon values ranging from [0.1, 0.18, 0.2, 0.22, 0.25, 0.3, 0.4] and submit the batch job using the Slurm script.
3. Finally, to gather the Rg values for different IDPs, run collect_Rg_data.py. This will generate subplots with all the Rg values for different epsilon.

If you have any questions or difficulty running, please contact:
Swarnadeep Seth: Swarnadeep.Seth@ucf.edu
Aniket Bhattacharya: Aniket.Bhattacharya@ucf.edu
