import os

ParticleN = 67
T_Kelvin=298
BOX_L = 600
ionic_concentration = 42 # in mM

os.system("python3 fasta2_three_letter.py") # create the three letter codes
epsilon_array = [0.1, 0.18, 0.2, 0.22, 0.25, 0.3, 0.4]

for i in epsilon_array:
    cmd = "mkdir -p " + "E_" + str(i)
    os.system (cmd)
    
    with open('code/ProteinBulk.py', 'r') as filee :
        filedata = filee.readlines()

        filedata[9] = f'ParticleN={ParticleN}' + '\n'
        filedata[12] = f'BOX_L = {BOX_L}' + '\n'
        filedata[17] = f'T_Kelvin={T_Kelvin}' + '\n'

        filedata[24] = f'EPSILON = {i} # KCal/mol' + '\n'
        filedata[26] = f'ionic_concentration = {ionic_concentration}*1e-3 # in M or mol/L' + '\n'

    with open(f'E_{i}/ProteinBulk.py', 'w') as filee:
        filee.writelines(filedata)

    filee.close()
    
    cmd=f'cp -a code/{ParticleN} ' + 'E_' + str(i)
    os.system(cmd)

    cmd='cp -a code/ThreeLetterCode.dat ' + 'E_' + str(i)
    os.system(cmd)

    cmd='cp -a code/stats_module.dat ' + 'E_' + str(i)
    os.system(cmd)

    cmd='cp -a code/chain_param.dat ' + 'E_' + str(i)
    os.system(cmd)

    cmd='cp -a code/long_config_analysis.py ' + 'E_' + str(i)
    os.system(cmd)
    
    cmd='cp -a code/run.sh ' + 'E_' + str(i)
    os.system(cmd)

    os.chdir(f"E_{i}")
    cmd='sbatch run.sh'
    os.system(cmd)

    os.chdir("../")