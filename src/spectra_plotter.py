import os, sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def make_df(LINES, cols):
    df_entries = {}
    for index, line in enumerate(LINES):
        linelist = line.split()
        df_entries[index] = linelist
    df = pd.DataFrame.from_dict(df_entries, orient="index", columns=cols)
    
    return df

def make_plot(df, ID):
    x = df["lambda"].apply(eval)
    y = df["energy"].apply(eval)                                                                                                                                                                           
    font = 'Arial'
    fig, ax = plt.subplots(figsize=(18, 9))                                                                                                                                                                
    plt.plot(x, y)
    plt.title("Absorption Spectra", fontsize=16, fontname=font)
    plt.xlabel(f"Wavelength ()", fontsize=12, fontname=font)
    plt.ylabel(f"Energy ()", fontsize=12, fontname=font)               
    fig.savefig("data/spectra_pngs/"+f"{ID}.png", format='png', dpi=450, bbox_inches='tight')
    plt.close(fig)

try:
    out = sys.argv[2] # supply as second arg if needed
except:
    out = "data/spectra_outputs/" # path to collection of tddft outputs
    
for spectra_file in os.listdir(out):

    with open(out+spectra_file, "r") as f:
        lines = f.readlines()
    for line_number, line_details in enumerate(lines):
        if "----" in line_details:
            if "ABSORPTION" in lines[line_number+1]:
                if "ELECTRIC" in lines[line_number+1]:
                    electric_DM = lines[line_number+6:line_number+26]
                    velocity_DM = lines[line_number+33:line_number+53]
                    CD = lines[line_number+60:line_number+80]
                    CD_velocity = lines[line_number+87:line_number+107]

    # make dataframes    
    df_electric_DM = make_df(electric_DM, ['state', 'energy', 'lambda', 'frequency_oscillation', 'T2', 'TX', 'TY', 'TZ'])
    df_velocity_DM = make_df(velocity_DM, ['state', 'energy', 'lambda', 'frequency_oscillation', 'P2', 'PX', 'PY', 'PZ'])
    df_CD = make_df(CD, ['state', 'energy', 'lambda', 'R', 'MX', 'MY', 'MZ'])
    df_CD_velocity = make_df(CD_velocity, ['state', 'energy', 'lambda', 'R', 'MX', 'MY', 'MZ'])

    if sys.argv[1] == "multiple":
        code = "_".join(spectra_file.split("_")[:3])
        dir = "mult_conf_dimer_spectra_csvs"
    elif sys.argv[1] == "single":
        code = "_".join(spectra_file.split("_")[:2]) + "_1"
        dir = "single_conf_dimer_spectra_csvs"
    os.system(f"mkdir -p {dir}")
    df_electric_DM.to_csv(f"{dir}/{code}_electric.csv")

'''
    try:
        if sys.argv[1] == 'monomers':
            ID = spectra_file.split("_")[0]            
            try:
                float(spectra_file.split("_")[1])
                df_electric_DM.to_csv("data/excited_state_monomer_spectra_csvs/"+f"{ID}_electric.csv")
                print(f"{ID} done")
            except:
                float(spectra_file.split("_")[1])
                df_electric_DM.to_csv("data/excited_state_monomer_spectra_csvs/"+f"{ID}_electric.csv")
                #df_velocity_DM.to_csv("data/spectra_csvs/"+f"{ID}_velocity.csv")
                #df_CD.to_csv("data/spectra_csvs/"+f"{ID}_CD.csv")
                #df_CD_velocity.to_csv("data/spectra_csvs/"+f"{ID}_CD_velocity.csv")
                #make_plot(df_electric_DM, ID)
                #make_plot(df_velocity_DM, ID)
                #make_plot(df_CD, ID)
                #make_plot(df_CD_velocity, ID)
                print(f"{ID} done")
        elif sys.argv[1] == 'dimers':
            ID = spectra_file.split("_")[0] + '_' + spectra_file.split("_")[1]        
            try:
                float(spectra_file.split("_")[1])
                df_electric_DM.to_csv("data/excited_state_dimer_spectra_csvs/"+f"{ID}_electric.csv")
                #df_velocity_DM.to_csv("data/spectra_csvs/"+f"{ID}_velocity.csv")
                #df_CD.to_csv("data/spectra_csvs/"+f"{ID}_CD.csv")
                #df_CD_velocity.to_csv("data/spectra_csvs/"+f"{ID}_CD_velocity.csv")
                #make_plot(df_electric_DM, ID)
                #make_plot(df_velocity_DM, ID)
                #make_plot(df_CD, ID)
                #make_plot(df_CD_velocity, ID)
                print(f"{ID} done")
            except:
                pass            
    except:
        ID = spectra_file.split("_")[0] + '_' + spectra_file.split("_")[1]
        try:
            float(spectra_file.split("_")[1])
            df_electric_DM.to_csv("data/spectra_csvs/"+f"{ID}_electric.csv")
            df_velocity_DM.to_csv("data/spectra_csvs/"+f"{ID}_velocity.csv")
            df_CD.to_csv("data/spectra_csvs/"+f"{ID}_CD.csv")
            df_CD_velocity.to_csv("data/spectra_csvs/"+f"{ID}_CD_velocity.csv")
            make_plot(df_electric_DM, ID)
            make_plot(df_velocity_DM, ID)
            make_plot(df_CD, ID)
            make_plot(df_CD_velocity, ID)
            print(f"{ID} done")
        except:
            pass      
'''  
