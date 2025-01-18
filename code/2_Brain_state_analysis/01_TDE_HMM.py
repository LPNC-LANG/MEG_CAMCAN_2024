
from sys import argv

if len(argv) != 3:
    print(
        "Please pass the number of states and run id, e.g. python 2_train_hmm.py 8 1"
    )
    exit()

n_states = int(argv[1])
run = int(argv[2])

import pickle
from glob import glob
from osl_dynamics.data import Data
from osl_dynamics.models.hmm import Config, Model

data = Data("/silenus/PROJECTS/pr-deepneuro/guichetc/MEG_data_TDE_HMM",
            use_tfrecord=True,
            buffer_size=1000,
            n_jobs=16)


config = Config(
    n_states=n_states,
    n_channels=data.n_channels,
    sequence_length=2000,
    learn_means=False,
    learn_covariances=True,
    batch_size=16,
    learning_rate=0.001,
    n_epochs=20
)

# Create model
model = Model(config)
model.summary()

## Training

# Initialisation
init_history = model.random_state_time_course_initialization(
    data, n_init=3, n_epochs=1
    )

from osl_dynamics import run_pipeline
# Full training
history = model.fit(data)

# Save trained model
model_dir = f"results/{n_states:02d}_states/run{run:02d}/model"
model.save(model_dir)

# Calculate the free energy
free_energy = model.free_energy(data)
history["free_energy"] = free_energy

## Save training history
with open(f"{model_dir}/history.pkl", "wb") as file:
    pickle.dump(history, file)

with open(f"{model_dir}/loss.dat", "w") as file:
    file.write(f"ll_loss = {history['loss'][-1]}\n")
    file.write(f"free_energy = {free_energy}\n")