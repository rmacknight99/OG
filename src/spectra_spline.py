import pandas as pd
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import daltonproject as dp 

parser = argparse.ArgumentParser(description='Fit Spectra Data')
parser.add_argument('--ID', type=str, 
                    help='ID to compare for ')
args = parser.parse_args()

if args.ID is None:
	raise Exception('No ID given')


# This script assumes
# Experimental data oliveira_gomes/data/experimental_spectras
# Computational data oliveira_gomes/data/excited_state_dimer_spectra_csvs
# We compare electric excited states 



# For some reason ... experimental data is named either is ID1-ID2 or ID1_ID2 so we test for either


ID = args.ID



experimental_path = 'data/experimental_spectras/' + args.ID + '.csv'
comp_path = 'data/excited_state_dimer_spectra_csvs/' + args.ID + '_electric.csv'


exp_data = pd.read_csv(experimental_path, names = ['wave length', 'absorption'])





data = pd.read_csv(comp_path)

x = data['lambda'].to_numpy()
y  = data['frequency_oscillation'].to_numpy()

print(x,y)
# x is a decreasing array
stop = np.argmin(x >= 300)
start = np.argmax((x <= 800) & (x > 0))
print(start, stop)
if start > stop:
	stop = len(x)

x = x[start:stop]
y = y[start:stop]

rows,  = x.shape
x = x.reshape((rows,1))
y = y.reshape((rows,1))

# The excited state calculation could have lambas between 100 - 2000 nm if not more .... we restrict the range to 300 - 800 nm to allow for some noise (if available)
# in the spline surrounding the visible region of interest (390 - 700)

#spl = UnivariateSpline(x, y)

#wave_lengths = np.arange(390, 701)

#absorption = spl(x)

#Normalize
#absorption /= sum(absorption)


# Wavelengths are in nm, nm -> eV is a 1239.8 

x = np.reciprocal(x)
x *= 1239.8

start = 1239.8 / 390
stop = 1239.8 / 700

e_range = np.linspace(start = start, stop = stop, num = 1000)

lambda_range = np.reciprocal(e_range)
lambda_range *= 1239.8
spectrum = dp.spectrum.convolute_one_photon_absorption(excitation_energies = x,
                                    oscillator_strengths = y,
				    energy_range = e_range)

print(spectrum)
#spectrum -= min(spectrum)
#spectrum /= (max(spectrum) - min(spectrum))



absorption = exp_data['absorption'].to_numpy()
absorption -= min(absorption)
absorption /= (max(absorption) - min(absorption))

plt.plot(lambda_range, spectrum, label = 'orca')
#plt.plot(exp_data['wave length'].to_numpy(), absorption, label = 'experimental')

plt.legend()
plt.show()

exit()




# Save plot in same format as orca_uv
#plt.scatter(x, absorption, label = 'spline')


plt.scatter(exp_data['wave length'].to_numpy(), exp_data['absorption'].to_numpy(), label = 'experimental')
plt.scatter(x, y, label = 'orca')
plt.title("Absorption Spectra", fontsize=16, fontname='Arial')
plt.xlabel(f"Wavelength ()", fontsize=12, fontname='Arial')
plt.ylabel(f"Predicted Absorption ()", fontsize=12, fontname='Arial')
plt.legend()
plt.savefig("test.png", format='png', dpi=450, bbox_inches='tight')
plt.close()



