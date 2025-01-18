"""Print free energy.

"""

# from sys import argv

# if len(argv) != 2:
#     print("Please pass the number of states, e.g. python 2_print_free_energy.py 8")
#     exit()
# n_states = int(argv[1])

import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

BASE_DIR = "E:/Research_Projects/MEG_CamCAN"

os.system("cls")

free_energy_data = {
    'state': [],
    'run': [],
    'free_energy': []
}

best_fe = np.Inf
best_state = None
best_run = None

for n_states in (6, 8):
    for run in range(1, 6):
        model_dir = f"{BASE_DIR}/TDE_HMM/results/{n_states:02d}_states/run{run:02d}/model"
        try:
            history = pickle.load(
                open(f"{model_dir}/history.pkl", "rb")
            )
            free_energy = history["free_energy"]
            print(f"run {run}: {free_energy}")

            free_energy_data['state'].append(n_states)
            free_energy_data['run'].append(run)
            free_energy_data['free_energy'].append(free_energy/1e5)

            if free_energy < best_fe:
                best_run = run
                best_fe = free_energy
                best_state = n_states
        except:
            print(f"run {run} missing")

print()
print("best state:", best_state, 
      "\nbest run:", best_run, 
      "\nwith free energy:", best_fe)

df = pd.DataFrame(free_energy_data)

# Create the violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x='state', y='free_energy', data=df, inner='quartile')
plt.title('Free Energy Distribution by State and Run')
plt.xlabel('Number of States')
plt.ylabel('Free Energy')
plt.show()