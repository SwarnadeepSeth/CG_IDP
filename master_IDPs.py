import os
import numpy as np
import pandas as pd

df  = pd.read_csv("idp_bank.csv")
df.columns = ["Name", "Length", "Cs", "Temp", "Sequence"]
print (df)

print ("="*60)
start = int(input("Enter the start index of the IDP to run: "))
end = int(input("Enter the end index of the IDP to run: "))

df = df.loc[start:end]
df = df.reset_index(drop=True)
print ("="*60)
print ("Creating the Run folders for the followings:")
print (df)
print ("="*60)

BOX_L = 600

for i in range (len(df)):
    idp_name = df["Name"][i]
    os.mkdir(idp_name)

    ParticleN = df['Length'][i]
    ionic_concentration = df['Cs'][i]
    T_Kelvin = df['Temp'][i]
    Seq = df["Sequence"][i]

    with open('master.py', 'r') as filee :
        filedata = filee.readlines()

        filedata[2] = f'ParticleN={ParticleN}' + '\n'
        filedata[3] = f'T_Kelvin={T_Kelvin}' + '\n'
        filedata[4] = f'BOX_L = {BOX_L}' + '\n'
        filedata[5] = f'ionic_concentration = {ionic_concentration} # in mM' + '\n'

    with open(f'{idp_name}/master.py', 'w') as filee:
        filee.writelines(filedata)

    os.system (f"cp -r code {idp_name}")
    os.system (f"cp fasta2_three_letter.py {idp_name}")

    fi = open(f"{idp_name}/{idp_name}.fasta", "w")    
    print (Seq, file=fi)
    fi.close()

