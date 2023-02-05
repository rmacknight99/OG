import pandas as pd
import argparse, os, itertools, sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import daltonproject as dp  
from sklearn.metrics import r2_score
import torch
import gpytorch
import moviepy.video.io.ImageSequenceClip
from rdkit.Chem import AllChem as Chem
import warnings
warnings.filterwarnings("ignore")
from morfeus_utils import *

# Identifying matches - for the computational data we have many conformers (we should Boltzmann weight the spectra?)
# For now this function will only operate on the first conformer

def identify_matches(single_path="comp_spectral_data/single_conf_dimer_spectra_csvs/", 
                     multiple_path="comp_spectral_data/mult_conf_dimer_spectra_csvs/", 
                     exp_path="exp_spectral_data/", conf_tracker={}):
    # paths
    single_conformers_path = single_path
    multiple_conformers_path = multiple_path
    experimental_path = exp_path

    # first find comp codes and files
    for f in os.listdir(single_conformers_path):
        conf_id = "1"
        code = sorted([int(i) for i in f.split("_")[:2]])
        x = "_".join(map(str, code))
        conf_tracker[x]["comp_file"].append(single_conformers_path+f)
    for f in os.listdir(multiple_conformers_path):
        conf_id = int(f.split("_")[2])
        code = sorted([int(i) for i in f.split("_")[:2]])
        x = "_".join(map(str, code))
        conf_tracker[x]["comp_file"].append(multiple_conformers_path+f)
    # now find the exp codes and files 
    sorted_exp_codes = []
    for f in os.listdir(exp_path):
        if f.endswith(".csv"):
            substring = f.split('.')[0]
            code = sorted([int(i) for i in substring.split('_')])
            x = "_".join(map(str, code))
            conf_tracker[x]["exp_file"] = exp_path+f  
    # now find the matches
    matches = []
    keys_to_remove = []
    for key, value in conf_tracker.items():
        if value["exp_file"] != "" and len(value["comp_file"]) > 0:
            matches.append(key)
        elif value["exp_file"] == "" or len(value["comp_file"]) == 0:
            keys_to_remove.append(key)
    #[conf_tracker.pop(key) for key in keys_to_remove]
    print(f"{len(matches)} matching comp and exp spectra")
    
    return matches, conf_tracker

# Lets first analze the differences between our experimental and computational spectras

# Get spectra from csv file for exp and comp

def get_data(data_path = None):
    if "exp" in data_path:
        data = pd.read_csv(data_path, names = ['wavelength', 'absorption'])
        x = data["wavelength"].tolist()[::-1]
        y = data["absorption"].tolist()
    else:
        data = pd.read_csv(data_path)
        x = data["lambda"].tolist()[::-1]
        y = data["frequency_oscillation"].tolist()
    y = [(h - min(y))/(max(y) - min(y)) for h in y][::-1] # min max scale y axis

    return x, y

# A function for getting a full computational spectra

def gen_spectrum_from_comp(x, y, n_pts = 311, start_nm = 390, stop_nm = 700):

    x, y = np.array(x), np.array(y)
    start = np.argmin(x >= start_nm)
    if x[-1] < stop_nm:
        stop = len(x)
    else:
        stop = np.argmax(x >= stop_nm)
    x = x[start:stop]
    y = y[start:stop]
    rows,  = x.shape
    x, y = x.reshape((rows,1)), y.reshape((rows,1))
    x = np.reciprocal(x) * 1239.8
    start_e, stop_e = 1239.8 / start_nm, 1239.8 / stop_nm
    e_range = np.linspace(start = start_e, stop = stop_e, num = n_pts)
    lambda_range = np.reciprocal(e_range) * 1239.8
    spectrum = dp.spectrum.convolute_one_photon_absorption(excitation_energies = x,
                                        oscillator_strengths = y,
                                        energy_range = e_range)
    spectrum -= min(spectrum)
    spectrum /= (max(spectrum) - min(spectrum))
    wavelength = np.linspace(start = start_nm, stop = stop_nm, num = n_pts)
    
    return wavelength, spectrum

# Boltzmann Weight Spectrum

def boltzmann_weight_spectrum(collection, weights=None):
    total_spectrum = np.zeros(311)
    for i, comp_file in enumerate(collection):
        comp_x, comp_y = get_data(data_path=comp_file)
        WL, spectrum = gen_spectrum_from_comp(comp_x, comp_y, n_pts=311)
        if weights is not None:
            try:
                total_spectrum += spectrum * weights[i]
            except:
                total_spectrum = total_spectrum
        else:
            total_spectrum += spectrum
    if weights is not None:
        FR = total_spectrum
        return WL, FR
    else:
        FR = total_spectrum / len(collection)
        return WL, FR

# A function to plot the spectra

def gen_plot(x, y, data_label, title, x_label="X", y_label="Y", color="black", marker=None, multiple=False, scatter=False):
    # make figure and axes
    fig, ax = plt.subplots(figsize=(12, 12))
    # aesthetic settings for plot
    #font = 'Arial'
    SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 16, 20, 24
    ax.grid(True, linewidth=1.0, color='0.95')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_axisbelow(True)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3.0)
    for tick in ax.get_yticklabels():
        #tick.set_fontname(font)
        tick.set_fontsize(SMALL_SIZE)
    for tick in ax.get_xticklabels():
        #tick.set_fontname(font)
        tick.set_fontsize(SMALL_SIZE)

    if multiple:
        for i in range(len(y)):
            plt.plot(x, y[i], label = data_label[i], color=color[i], marker=marker, linewidth=6.0)
    elif scatter:
        plt.scatter(x, y, label=data_label, color=color, marker=".", linewidth=6.0)
    else:
        plt.plot(x, y, label=data_label, color=color, marker=marker, linewidth=6.0)
    ax.legend(fontsize=SMALL_SIZE, loc=0)
    plt.title(title, fontsize=BIGGER_SIZE)
    plt.xlabel(f"{x_label}", fontsize=MEDIUM_SIZE)
    plt.ylabel(f"{y_label}", fontsize=MEDIUM_SIZE)
    plt.xticks([400, 500, 600, 700])
    plt.close()
    return fig, ax

# A function for getting some stats when comparing comp and exp spectra

def get_stats(exp_y, comp_y):

    exp_y = np.asarray(exp_y)
    comp_y = np.asarray(comp_y)
    abs_diff = np.sum(np.abs(exp_y - comp_y))
    n = exp_y.shape[0]

    MAE = abs_diff/(n)
    r2 = r2_score(exp_y, comp_y)
    if r2 > 1 or r2 < 0:
        r2 = 0
        
    return MAE, r2

# A function for comparing the spectra of matches

def compare(conf_tracker={}, boltz_weight=False):

    keys_to_remove = []
    data_dict = {"comp_file": [], "exp_file": [], "MAE_before": [], "R2_before": [], "WL": [], "comp_y": [], "exp_y": []}
    for complex, information in conf_tracker.items():
        print(f"\tWorking on {complex}...")
        try:
            comp_file = information["comp_file"]
        except:
            pass
        if comp_file != []:
            data_dict["comp_file"].append(comp_file)
            if len(comp_file) > 1:
                if boltz_weight:
                    print("\t\tBoltzmann weighting the spectral data spectrum ...")
                    ensemble = load_ensemble('_'.join(comp_file[0].split("/")[-1].split("_")[:2])+"_")
                    _ = gen_crest(ensemble)
                    weights = get_boltz_weights(".")
                    print(f"\t\t\tweights = {', '.join([str(round(w, 2)) for w in weights])}")
                    WL, FR = boltzmann_weight_spectrum(comp_file, weights=weights)
                else:
                    WL, FR = boltzmann_weight_spectrum(comp_file)
            else:
                comp_x, comp_y = get_data(data_path=comp_file[0])
                WL, FR = gen_spectrum_from_comp(comp_x, comp_y, n_pts = 311)
            data_dict["WL"].append([round(i, 2) for i in WL])
            data_dict["comp_y"].append([round(i, 4) for i in FR])
        else:
            WL, FR = None, None
            data_dict["comp_file"].append(0)
            data_dict["WL"].append([])
            data_dict["comp_y"].append([])
        try:
            exp_file = information["exp_file"]
        except:
            pass
        if exp_file != "":
            data_dict["exp_file"].append(exp_file)
            exp_x, exp_y = get_data(data_path=exp_file)
            data_dict["exp_y"].append([round(i, 4) for i in exp_y])
            if FR is not None:
                MAE, r2 = get_stats(exp_y, FR)
                data_dict["MAE_before"].append(MAE)
                data_dict["R2_before"].append(r2)
                match_code = complex.split("_")
                fig, ax = gen_plot(WL, [exp_y, FR], ["experimental", "computational"], 
                                   f"Interpolated Comp Spectrum vs. Exp Spectrum\n monomers = {match_code}", 
                                   x_label="Wavelength (nm)", y_label="Normalized Intensity", color=["blue", "red"], multiple=True)
                os.system("mkdir -p before_training/")
                fig.savefig(f"before_training/comparison_{match_code[0]}_{match_code[1]}.png")
            else:
                data_dict["MAE_before"].append(None)
                data_dict["R2_before"].append(None)
        else:
            data_dict["MAE_before"].append(None)
            data_dict["R2_before"].append(None)
            data_dict["exp_file"].append(None)
            data_dict["exp_y"].append([])
    df = pd.DataFrame.from_dict(data_dict)
    df = df.drop(df[df.comp_file == 0].index)
    return df

# A function for loading monomer spectra

def load_monomer(monomer_path, n_pts=32):

    data = {}
    for csv in os.listdir(monomer_path):
        comp_file = monomer_path+csv
        comp_x, comp_y = get_data(data_path=comp_file)
        WL, FR = gen_spectrum_from_comp(comp_x, comp_y, n_pts = n_pts)
        data[comp_file] = [WL, FR]
    return data
    
# Gather a subset of 'similar' spectra

def get_similar_spectra(threshold, df):

    return df.where(df["R2_before"] > threshold).dropna()

# Training functions

# A function for splitting our data into train and test

def split_data(samples_df, test_size = 0.2):

    amt_of_data = samples_df.shape[0]

    train_indices = np.sort(np.random.choice(amt_of_data, int(amt_of_data*(1-test_size)), replace=False)) # generate indices for training data
    test_indices = np.asarray([i for i in range(amt_of_data) if i not in train_indices])

    training_df = samples_df.iloc[train_indices, :]
    testing_df = samples_df.iloc[test_indices, :]

    return training_df, testing_df

# Exact GP model with the RBF Kernel (single task, could predict a intensity given a wavelength for a complex)

class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel())

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

# Exact GP model with RBF Kernel (multi-task, could predict entire spectrum for a complex)

class MultitaskGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(MultitaskGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.MultitaskMean(gpytorch.means.ConstantMean(), num_tasks=list(train_y.size())[1])
        self.covar_module = gpytorch.kernels.MultitaskKernel(gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel()), num_tasks=list(train_y.size())[1], rank=1)

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultitaskMultivariateNormal(mean_x, covar_x)
        
# A function for generating the features we can use for training

def gen_features(df, monomer_path, monomer_data, monomer_features=False, comp_features=False, add_fps=False, fp_bits=64):

    features_col = []
    for index, row in df.iterrows():
        feats = []

        if comp_features:
            d_features = row["comp_y"].tolist()
            feats.extend(d_features)

        if monomer_features:
            cf = row["comp_file"][0]
            m1_file = monomer_path+cf.split("/")[-1].split("_")[0]+"_electric.csv"
            m2_file = monomer_path+cf.split("/")[-1].split("_")[1]+"_electric.csv"
            m1_features = monomer_data[m1_file][1].tolist() # just intensity
            m2_features = monomer_data[m2_file][1].tolist() # just intensity
            feats.extend(m1_features)
            feats.extend(m2_features)


        if add_fps:
            smiles_df = pd.read_csv("SMILES.csv", index_col=0)
            smiles_lookup = smiles_df.to_dict(orient="index")
            cf = row["comp_file"][0]
            monomers = cf.split("/")[-1].split("_")
            m1, m2 = int(monomers[0]), int(monomers[1])
            s1, s2 = smiles_lookup[m1]["SMILES"], smiles_lookup[m2]["SMILES"]
            fp1 = list(Chem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s1), 2, nBits=fp_bits))
            fp2 = list(Chem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s2), 2, nBits=fp_bits))
            feats.extend(fp1)
            feats.extend(fp2)

        features_col.append(feats)
    try:
        df.insert(2, "training_features", features_col)
    except ValueError:
        df["training_features"] = features_col

    return df

# A function for training a multitask GP

def train_multitask_GP(training_df, training_data="training_features", target="exp_y", n_iterations = 250, GPU=False, train=True):

    train_x = torch.tensor(training_df[training_data].tolist()).float()
    train_y = torch.tensor(training_df[target].tolist()).float()

    likelihood = gpytorch.likelihoods.MultitaskGaussianLikelihood(num_tasks = train_y.size()[1]) # initialize the Multitask Gaussian likelihood with as many tasks as elements in each output
    model = MultitaskGPModel(train_x, train_y, likelihood) # intialize Multitask GP model 

    if GPU:
        print("**Using GPU**")
        train_x = train_x.cuda()
        train_y = train_y.cuda()
        model = model.cuda()
        likelihood = likelihood.cuda()

    model.train()
    likelihood.train()

    optimizer = torch.optim.Adam(model.parameters(), lr=0.1) # adam optimizer
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model) # Marginal Log Likelihood loss
    if train:
        print("\n----- Training Model -----\n")
        for i in range(n_iterations): # for a number of training iterations
            optimizer.zero_grad() # zero all tensor gradients
            output = model(train_x) # compute output of model on training data
            loss = -mll(output, train_y) # compute loss
            loss.backward() # back prop pgradients
            if (i+1) % 50 == 0 or (i+1) == n_iterations:
                print(f"\tTraining iteration {i+1}/{n_iterations} ---> "+" loss: {:.3f}".format(loss)) 
            optimizer.step() # step with optimizer
        print("\nFinal Loss: {:.3f}".format(loss))
    else:
        print("\n----- Initializing and returning model and likelihood -----\n")

    model.eval() # turn model into evaluation mode
    likelihood.eval() # turn likeihood into evaluation mode

    return model, likelihood # return model and likelihood in evaluation mode

# A function for predicting from a multitask GP model

def predict(model, likelihood, testing_df, 
            training_data="training_features", target="exp_y", 
            train=False, compare=False, plot=False, interval=20, GPU=False, load=None):
    
    if load is not None:
        state_dict = torch.load(load)
        model.load_state_dict(state_dict)
    if GPU:
        print("**Using GPU**")
    x_data = testing_df[training_data].tolist()
    y_data = testing_df[target].tolist()  
    mean_preds = []
    print("\n----- Running Predictions -----\n")
    if plot:
        WL = testing_df["WL"].tolist()[0]
        files = [i[0].split("/")[-1] for i in testing_df["comp_file"].tolist()]
        match_codes = [[j.split("_")[0], j.split("_")[1]] for j in files]

    if compare:
        r2_scores = []
        MAEs = []

    for i, x in enumerate(x_data):
        if i == 0 or i % interval == 0:
            print(f"\tlabeling {i+1}/{testing_df.shape[0]}, train: {train}")
        elif (i+1) == len(x_data):
            print(f"\tlabeling {i+1}/{testing_df.shape[0]}, train: {train}")
        x = torch.tensor([x]).float()
        if GPU:
            x = x.cuda()
        pred = likelihood(model(x)) # get hyperparameter predictions
        mean_pred = pred.mean.tolist() # get mean predictions
        mean_preds.append([round(i, 3) for i in mean_pred[0]])
        if compare:
            MAE, r2 = get_stats(mean_pred[0], y_data[i])
            r2_scores.append(r2)
            MAEs.append(MAE)

        if plot:

            if train:
                fig, ax = gen_plot(WL, [y_data[i], mean_pred[0]], ["experimental", "computational"], 
                            f"Interpolated Comp Spectrum vs. Exp Spectrum\n monomers = {match_codes[i]}, train", 
                            x_label="Wavelength (nm)", y_label="Normalized Intensity", color=["blue", "red"], multiple=True)
                os.system("mkdir -p after_training/")
                fig.savefig(f"after_training/comparison_{match_codes[i][0]}_{match_codes[i][1]}_train.png")

            else:
                fig, ax = gen_plot(WL, [y_data[i], mean_pred[0]], ["experimental", "predicted"], 
                            f"Predicted Spectrum vs. Exp Spectrum\n monomers = {match_codes[i]}, test", 
                            x_label="Wavelength (nm)", y_label="Normalized Intensity", color=["blue", "green"], multiple=True)
                os.system("mkdir -p after_training/")
                fig.savefig(f"after_training/comparison_{match_codes[i][0]}_{match_codes[i][1]}_test.png")

    if compare:
        testing_df["R2_after"] = r2_scores
        testing_df["MAE_after"] = MAEs
        r2_diff = np.mean(testing_df["R2_after"] - testing_df["R2_before"])
        mae_diff = np.mean(testing_df["MAE_before"] - testing_df["MAE_after"])
        print(f"Average increase in R2 = {round(r2_diff, 4)}, Average decrease in MAE = {round(mae_diff, 4)}")

    return mean_preds, testing_df

def make_gif(images, gif_name="my.gif", length=45):

    # make clip with desired length in seconds and save 
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(list(images), fps=len(images)//length)
    clip.write_gif(gif_name)

def gather_images(directory):
    images = []
    for i in os.listdir(directory):
        if i.endswith(".png"):
            images.append(directory+i)
    return images

def shrink_spectrum(df, cut=20, loc="comp_y"):
    
    vals = df[loc].tolist()
    new_vals = []
    for data in vals:
        new_data = []
        for i, e in enumerate(data):
            if i % cut == 0:
                new_data.append(e)
        new_vals.append(new_data)
    
    df[loc] = new_vals
    return df