from model_utils import *
from morfeus_utils import *

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
                        data, monomer_features=True, add_fps=True, fp_bits=196)
testing = gen_features(testing, "comp_spectral_data/excited_state_monomer_spectra_csvs/",
                       data, monomer_features=True, add_fps=True, fp_bits=196)

feature_size = len(training["training_features"].tolist()[0])
training, testing = shrink_spectrum(training, cut=10), shrink_spectrum(testing, cut=10)
target_size = len(training["comp_y"].tolist()[0])

print(f"\n  training examples: {training.shape[0]}\n  testing examples: {testing.shape[0]}\n  feature size: {feature_size}\n  target size: {target_size}\n")

# train a multitask GP to predict COMPUTATIONAL DIMER SPECTRUM from features
comp_model, comp_likelihood = train_multitask_GP(training, target="comp_y", n_iterations=250, GPU=True)
os.system("rm -rf after_training/")

# Run prediction and replace COMP spectra with prediction from MODEL
mean_train_preds, training = predict(comp_model, comp_likelihood, training, target="comp_y", train=True, interval=100)
mean_test_preds, testing = predict(comp_model, comp_likelihood, testing, target="comp_y", interval=50)

training["comp_y"] = mean_train_preds
testing["comp_y"] = mean_test_preds
