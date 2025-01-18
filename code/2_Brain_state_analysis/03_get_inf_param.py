"""Get inferred parameters.

"""

import os
# from sys import argv

# if len(argv) != 3:
#     print("Please pass the number of states and run id, e.g. python 4_get_inf_params.py 8 1")
#     exit()

# n_states = int(argv[1])
# run = int(argv[2])


n_states = 8
run = 4

BASE_DIR = "/media/clement/Transcend_Clement2/Research_Projects/MEG_CamCAN"

model_dir = f"{BASE_DIR}/TDE_HMM/results/{n_states:02d}_states/run{run:02d}/model"
output_dir = f"{BASE_DIR}/TDE_HMM/results/inf_params/{n_states:02d}_states"

os.makedirs(output_dir, exist_ok=True)

print("Importing packages")
import pickle
import numpy as np
from osl_dynamics.data import Data
from osl_dynamics.models import load

# Load training data
dataset = Data(
    f"{BASE_DIR}/TDE_HMM/TDE_data",
    n_jobs=16,
    use_tfrecord=True
)

# Load model
model = load(model_dir)
model.summary()

# Get inferred alphas (state probabilities)
alpha = model.get_alpha(dataset)
pickle.dump(alpha, open(f"{output_dir}/alp.pkl", "wb"))

# Get inferred covariances
covs = model.get_covariances()
np.save(f"{output_dir}/covs.npy", covs)