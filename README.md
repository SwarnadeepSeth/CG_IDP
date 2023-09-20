# CG_IDP (Coarse-grained model of IDPs)

The code uses a hydropathy scale called M3 for each of the amino acids, which was developed by Tesei et al. in 2021 ([DOI Link](https://doi.org/10.1073/pnas.2111696118){:target="_blank"}). 
To run the code for 33 IDPS (listed in idp_bank.csv), please follow the following procedure:
1. Run master_IDPs.py, and when prompted, specify the index of the IDP sequence to run. This will create a directory with the name of the IDP.
2. cd to the directory and run master.py. This will generate folders with epsilon values ranging from [0.1, 0.18, 0.2, 0.22, 0.25, 0.3, 0.4] and submit the batch job using the Slurm script.
3. Finally, to gather the Rg values for different IDPs, run collect_Rg_data.py. This will generate subplots with all the Rg values for different epsilons.

If you have any questions or difficulty running, please contact: <br>
Swarnadeep Seth: Swarnadeep.Seth@ucf.edu and <br>
Aniket Bhattacharya: Aniket.Bhattacharya@ucf.edu
