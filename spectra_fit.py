from model_utils import *
from morfeus_utils import *

torch.set_num_threads(int(os.cpu_count()))
# initialize the conformer tracker
monomers = np.arange(50)+1
combs = []
for comb in itertools.combinations(monomers, 2):
    combs.append(list(comb))
conf_tracker = dict(("_".join(map(str, el)), {"exp_file": "", "comp_file": []}) for el in combs)

# identify matches with experimental data
matches, conf_tracker = identify_matches(conf_tracker=conf_tracker)
avg_confs = np.mean(np.array([len(conf_tracker[key]) for key in list(conf_tracker.keys())]))
print(f"on average each complex has {round(avg_confs, 3)} conformers")

##### PREDICTING DIMER COMPUTATIONAL SPECTRUM FROM MONOMER SPECTRA AND FINGERPRINTS #####

# First load all complexes
master = pd.read_pickle("complexes_1168_boltz_weighted.pkl")
samples = master

# Load monomer data
data = load_monomer("comp_spectral_data/excited_state_monomer_spectra_csvs/", n_pts=64)

# Split into train and test sets 
training, testing = split_data(samples, test_size=0.25)

# Generate features from monomer spectra and fingerprints for both train and testing sets
training = gen_features(training, "comp_spectral_data/excited_state_monomer_spectra_csvs/", 
                        data, monomer_features=True, add_fps=True, fp_bits=256)
testing = gen_features(testing, "comp_spectral_data/excited_state_monomer_spectra_csvs/",
                       data, monomer_features=True, add_fps=True, fp_bits=256)

feature_size = len(training["training_features"].tolist()[0])
training, testing = shrink_spectrum(training, cut=1), shrink_spectrum(testing, cut=1)
target_size = len(training["comp_y"].tolist()[0])

print(f"\n  training examples: {training.shape[0]}\n  testing examples: {testing.shape[0]}\n  feature size: {feature_size}\n  target size: {target_size}")

# train a multitask GP to predict COMPUTATIONAL DIMER SPECTRUM from features
GPU = torch.cuda.is_available()
if not os.path.exists("comp_model.pth"):
    train = True
else:
    train = False
comp_model, comp_likelihood = train_multitask_GP(training, target="comp_y", n_iterations=300, GPU=GPU, train=train)
os.system("rm -rf after_training/")
if train:
    torch.save(comp_model.state_dict(), "comp_model.pth")

# Run prediction and replace COMP spectra with prediction from MODEL
comp_model 
mean_train_preds, training = predict(comp_model, comp_likelihood, training, target="comp_y", train=True, interval=100, GPU=GPU, load="comp_model.pth")
mean_test_preds, testing = predict(comp_model, comp_likelihood, testing, target="comp_y", interval=50, GPU=GPU, load="comp_model.pth")

training["comp_y"] = mean_train_preds
testing["comp_y"] = mean_test_preds


# Rejoin training and testing DataFrames 
exp_samples = pd.concat([training, testing]).dropna()

# Split into train and test sets 
training, testing = split_data(exp_samples, test_size=0.25)

# Generate features from dimer spectra, monomer spectra, and fingerprints for both train and testing sets
training = gen_features(training, "comp_spectral_data/excited_state_monomer_spectra_csvs/", 
                        data, monomer_features=True, add_fps=True, fp_bits=128, comp_features=True)
testing = gen_features(testing, "comp_spectral_data/excited_state_monomer_spectra_csvs/",
                       data, monomer_features=True, add_fps=True, fp_bits=128, comp_features=True)
                       
feature_size = len(training["training_features"].tolist()[0])
target_size = len(training["exp_y"].tolist()[0])

print(f"training examples: {training.shape[0]}\ntesting examples: {testing.shape[0]}\nfeature size: {feature_size}\ntarget size: {target_size}")

##### predicting EXPERIMENTAL DIMER SPECTRUM from PREDICTED COMP SPECTRA, MONOMER SPECTRA, and FPs #####
GPU = torch.cuda.is_available()
if not os.path.exists("exp_model.pth"):
    train = True
else:
    train = False
exp_model, exp_likelihood = train_multitask_GP(training, target="exp_y", n_iterations=250, GPU=GPU, train=train)
if train:
     torch.save(comp_model.state_dict(), "exp_model.pth")

# Run predictions
mean_train_preds, training = predict(exp_model, exp_likelihood, training, target="exp_y", train=True, plot=True, compare=True, load="exp_model.pth")
mean_test_preds, testing = predict(exp_model, exp_likelihood, testing, target="exp_y", plot=True, compare=True, load="exp_model.pth")
