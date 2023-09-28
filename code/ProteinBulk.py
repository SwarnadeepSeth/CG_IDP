import os
import hoomd
import hoomd.md
import gsd.hoomd
import numpy as np

hoomd.context.initialize("");

# ========================= System Parameters =======================================
ParticleN=71
dimen = 3
# Box Dimentions
BOX_L = 300

#kappa = 0.0

# ========================= Simulation Constants =======================================
T_Kelvin=298
# In kJ units
#KT = T_Kelvin*8.3145e-3 # Boltzmann constant in kJ/mol/K
#EPSILON = 0.2 * 4.184 # kJ/mol

# In kCal units
KT = T_Kelvin*0.001987204259 # Boltzmann constant in kCal/mol/K
EPSILON = 0.2 # KCal/mol

ionic_concentration = 199*1e-3 # in M or mol/L

fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3 #temperature dependent dielectric constant of water
epsw = fepsw(T_Kelvin) # dielectric constant of water at T 
lB = (1.6021766**2/(4*np.pi*8.854188*epsw))*(6.022*1000/KT)/4.182 # Bjerrum length in nm
#lB = (1.6021766**2/(4*np.pi*8.854188*80))*(6.022*1000/KT)/4.182 #  Bjerrum length in nm

yukawa_eps = lB*KT
yukawa_kappa = np.sqrt(8*np.pi*lB*ionic_concentration*6.022/10)

# ========================= System Parameters =======================================
fi_sys = open("System_Parameters.dat", "w")
fi_sys.write("Parameters and Constants ============================\n")
fi_sys.write("Number of particles: %d\n" % ParticleN)
fi_sys.write("Box length (nm): %f\n" % BOX_L)
fi_sys.write("Temperature (K): %d\n" % T_Kelvin)
fi_sys.write("LJ epsilon (kCal/mol): %f\n" % EPSILON)
fi_sys.write("Boltzmann constant*T (kCal/mol): %f\n" % KT)
fi_sys.write(f"Dielectric constant of water at T={T_Kelvin}: %f\n" % epsw)
fi_sys.write("Bjerrum length (nm): %f\n" % lB)
fi_sys.write("Ionic concentration (M): %f\n" % ionic_concentration)
fi_sys.write("Yukawa epsilon (kCal/mol): %f\n" % yukawa_eps)
fi_sys.write("Yukawa kappa (nm^-1): %f\n" % yukawa_kappa)
fi_sys.write("=================================================== \n")
fi_sys.close()
# ====================================================================================

TIME_STEP = 0.01
nwrite = 1000 # Frequency of saving configurations

Rouse_Steps = int((5*(ParticleN)**2.2)/TIME_STEP)
Run_Steps = Rouse_Steps
#Rouse_Steps = 1000000
#Run_Steps = 20000000

seed = np.random.randint(5000)
print ("Seed for the Run:", seed)

# ========================= Particles Connection Initialization =======================================
bond_gr=[]
for i in range(ParticleN-1):
	a_bond=[]
	a_bond.append(i)
	a_bond.append(i+1)
	bond_gr.append(a_bond)

angle_gr=[]	
for i in range(ParticleN-2):
	a_angle=[]
	a_angle.append(i)
	a_angle.append(i+1)
	a_angle.append(i+2)
	angle_gr.append(a_angle)
	
# ========================= System Initialization =======================================
system=gsd.hoomd.Snapshot()

aakeys = np.loadtxt('stats_module.dat', dtype=np.str, usecols=(0), unpack=True)
aakeys_mass = np.loadtxt('stats_module.dat', usecols=(1), unpack=True)
aakeys_chgs = np.loadtxt('stats_module.dat', usecols=(2), unpack=True)
aakeys_sigmas = np.loadtxt('stats_module.dat', usecols=(3), unpack=True)
aakeys_sigmas_scaled = aakeys_sigmas/10. # Angstrom to nm
aakeys_lambdas = np.loadtxt('stats_module.dat', usecols=(4), unpack=True)
aavalues = np.loadtxt('chain_param.dat', usecols=(1,2,3,4,5), unpack=True)

# Box description
system.configuration.dimensions = dimen
system.configuration.box = [BOX_L,BOX_L,BOX_L,0,0,0]

# Particle description
system.particles.N = ParticleN
system.particles.position = np.loadtxt(str(ParticleN))
system.particles.types = aakeys
system.particles.typeid = aavalues[0].astype(int)
system.particles.mass = aavalues[1]
system.particles.charge = aavalues[2]

# Bond description
system.bonds.N = ParticleN-1
system.bonds.types = ['AA_bond']
system.bonds.group = bond_gr

# Angle description
system.angles.N = ParticleN-2
system.angles.types = ['AA_angle']
system.angles.group = angle_gr

# Write intial chain config into gsd file
fi_start = gsd.hoomd.open(name='Start_Config.gsd', mode='wb')
fi_start.append(system)
#fi_start.close()

system = hoomd.init.read_gsd('Start_Config.gsd')

# ========================= Neighbor List =======================================
#nl = hoomd.md.nlist.cell(); # Change to tree for large box size
nl = hoomd.md.nlist.tree();

# ========================= Force Fields =======================================
# Custom Ashbaugh Potential
lj1 = hoomd.md.pair.lj(r_cut=3.5, nlist=nl, name="1")
lj2 = hoomd.md.pair.lj(r_cut=3.5, nlist=nl, name="2")

# Harmonic Bonds
harmonic=hoomd.md.bond.harmonic()
#harmonic.bond_coeff.set('AA_bond', k=8033, r0=0.38) # k in kJ/mol/nm^2 and r0 in nm
harmonic.bond_coeff.set('AA_bond', k=1920, r0=0.38) # k in kCal/mol/nm^2 and r0 in nm


# FENE Potential
#fene = hoomd.md.bond.fene()
#fene.bond_coeff.set('AA_bond', k=30.0, r0=1.5, sigma=1.0, epsilon= 0.0)
		
# Angle Potential: 
# V = k[1-cos(theta - theta0)], where theta0 = np.pi
# Force: T = - dV/d(theta)

#def bend_pot(theta, kappa):
#	V = kappa * (1.0+np.cos(theta));
#	T = kappa*np.sin(theta);
#	return (V,T)

#btable = hoomd.md.angle.table(width=1000)
#btable.angle_coeff.set('AA_angle', func=bend_pot, coeff=dict(kappa=kappa))

# Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut=3.5, nlist=nl)

# ========================= Pair Coefficients =======================================
for i in range (len(aakeys)):
	for j in range (len(aakeys)):
		sig_eff = 0.5*(aakeys_sigmas_scaled[i] + aakeys_sigmas_scaled[j]) 
		lambda_eff = 0.5*(aakeys_lambdas[i] + aakeys_lambdas[j]) 
		#print (i, j, aakeys[i], aakeys[j], sig_eff, lambda_eff)

		lj1.pair_coeff.set(aakeys[i], aakeys[j], epsilon=EPSILON*(1-lambda_eff), sigma=sig_eff, r_cut=(2**(1./6))*sig_eff)
		lj2.pair_coeff.set(aakeys[i], aakeys[j], epsilon=EPSILON*lambda_eff, sigma=sig_eff, r_cut=3.5)
  
		yukawa.pair_coeff.set(aakeys[i], aakeys[j], epsilon=yukawa_eps*aakeys_chgs[i]*aakeys_chgs[j], kappa=yukawa_kappa, r_cut=3.5)

lj1.set_params(mode='shift')
yukawa.set_params(mode='shift')
# nl.reset_exclusions(exclusions = []);

# ========================= MD Integrator =======================================
hoomd.md.integrate.mode_standard(dt=TIME_STEP);
all = hoomd.group.all();
integrator = hoomd.md.integrate.langevin(group=all, kT=KT, seed=seed, noiseless_t=False, noiseless_r=False)

for i,aa in enumerate(aakeys):
	GAMMA = aakeys_mass[i]/1000.
	integrator.set_gamma_r(aa, gamma_r=GAMMA)
	
# ========================= Warmup Run =======================================
hoomd.run(Rouse_Steps, quiet=True);
hoomd.dump.gsd("Equilibration_Config.gsd", period=nwrite, group=all, overwrite=True);
integrator.disable()

# ========================= MD Integrator =======================================
hoomd.md.integrate.mode_standard(dt=TIME_STEP);
all = hoomd.group.all();
integrator = hoomd.md.integrate.langevin(group=all, kT=KT, seed=seed, noiseless_t=False, noiseless_r=False)
for i,aa in enumerate(aakeys):
	GAMMA = aakeys_mass[i]/1000.
	integrator.set_gamma_r(aa, gamma_r=GAMMA)

# ========================= Print Thrmodynamic Quantities =======================================
#hoomd.analyze.log(filename="mddata.dat",quantities=['potential_energy','kinetic_energy','pair_lj_energy','bond_fene_energy', 'temperature'],period=10*nwrite,overwrite=True);
#hoomd.analyze.log(filename="mddata.dat",quantities=['potential_energy','kinetic_energy','pair_ashbaugh_energy','bond_harmonic_energy','pair_yukawa_energy', 'temperature'],period=10*nwrite,overwrite=True);
hoomd.analyze.log(filename="mddata.dat",quantities=['potential_energy','kinetic_energy','bond_harmonic_energy','pair_yukawa_energy', 'temperature'],period=10*nwrite,overwrite=True);

# ========================= Trajectory Print for Movie =======================================
hoomd.dump.gsd("Running_Config.gsd", period=nwrite, group=all, overwrite=True);
#hoomd.dump.dcd(filename="Running_Config.dcd", period=nwrite)

# ========================= Run Steps =======================================
hoomd.run(Run_Steps);

# ========================= Run Analysis =======================================
#os.system("python3 config_analysis.py")











 
