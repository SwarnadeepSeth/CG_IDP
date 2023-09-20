import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as scs

def categories(series):
    return range(int(series.min()), int(series.max()) + 1)

# ========================================================================================= #
# load the idp names
df= pd.read_csv('idp_bank.csv')
df.columns = ['Name', 'Length', 'Cs', 'Temp', 'Sequence']

# ========================================================================================= #
# Experiential data
Rg_exp = {'CspTm':1.47, 'IN':2.16, 'ProTaN':3.78, 'ProTaC':3.70, 'R15':2.33,
          'R17':2.29, 'hCyp':2.51, 'L':2.37, 'ACTR':2.51, 'hNHE1cdt':3.63,
          'sNase':2.12, 'synuclein':3.30, 'SIC1':3.21, 'HST5':1.38, 'CoINT':2.83,
          'p15PAF':2.81, 'FhuA':3.34, 'OPN':5.13, 'An16':4.40, 'p53':2.87, 
          'Nucleoporin153':2.4, 'ERMTADn': 3.96, 'SH4UD':2.82, 'K19':3.5, 'K18':3.8,
          'K17':3.6, 'K10':4.0, 'K27':3.7, 'K16':3.9, 'K25':4.1, 'K32':4.15,
          'K23':4.9, 'K44':5.2}

df['Rg_expt'] = df['Name'].map(Rg_exp)

# ========================================================================================= #
# Assign the ID
df_id = pd.read_csv('IDP_ID.csv', header=None)
df_id.columns = ['Name', 'ID']

# Assign the ID to df dataframe by comparing the Name column
df = pd.merge(df, df_id, on='Name')

# ========================================================================================= #
Epsilon_list = [0.1, 0.18, 0.2, 0.25, 0.3, 0.4]

col_num = 7
for eps in Epsilon_list:
    eps_array = []
    for i in df['Name']:

        filepath = i + f'/E_{eps}/config_stat.dat'

        # check if the file exists
        if not os.path.isfile(filepath):
            #print ("File path {} does not exist!".format(filepath))
            eps_array.append(np.nan)
        else:
            data = np.loadtxt(filepath, skiprows=1)
            eps_array.append(data[col_num])
            #print (i, eps, data[8])

    df[f'E_{eps}'] = np.array(eps_array)

# write the dataframe 
print (df)

# write the df into a csv
df.to_csv('all_Rg_data.csv', index=False)

# ========================================================================================= #
# plot the data into 5 column x 1 rows subplot for each epsilon
fig, ax = plt.subplots(1, 6, figsize=(16, 3), sharey=True)
for i, eps in enumerate(Epsilon_list):
    ax[i].plot(df['Rg_expt'], df[f'E_{eps}'], 'o', markersize=5, color='blue', alpha=0.5)
    # label the data points
    for j, txt in enumerate(df['ID']):
        ax[i].annotate(txt, (df['Rg_expt'][j], df[f'E_{eps}'][j]), fontsize=12)
    ax[i].set_title(fr'$\epsilon$ = {eps}', fontsize=20)

    ax[i].tick_params(axis='both', which='major', labelsize=15)
    ax[i].tick_params(axis='both', which='minor', labelsize=15)
    ax[i].set_xlim(1, 6.2)
    ax[i].set_ylim(1, 6.2)

    # plot 45 st degree line
    ax[i].plot([1, 6], [1, 6], 'k--', lw=2)
    #ax[i].grid(True)

    # show correlation coefficient at title
    #corr = df['Rg_expt'].corr(df[f'E_{eps}'])
    #ax[i].set_title(fr'$\epsilon$ = {eps}, corr = {corr:.2f}', fontsize=20)

    # calculate the MSE
    MSE = np.mean((df['Rg_expt'] - df[f'E_{eps}'])**2)
    legend1 = f'MSE = {MSE:.2f}'

    # calculate the chi2
    chi2 = np.mean((df['Rg_expt'] - df[f'E_{eps}'])**2 / df['Rg_expt']**2)

    legend2 = r'$\chi^2$ = ' + f'{chi2:.2f}'
    combined_legend = legend1 + '\n' + legend2
    ax[i].legend(title=legend1, title_fontsize=14, frameon=False)


#ax[2].set_xlabel(r'$R_g^{expt}$', fontsize=25)
ax[0].set_ylabel(r'$R_g^{M3}$', fontsize=25, color='blue')
plt.tight_layout()
plt.savefig('Rg_expt_vs_sim.pdf')
plt.show()


