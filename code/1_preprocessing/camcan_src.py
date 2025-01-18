"""Source reconstruction: forward modelling, beamforming and parcellation.

"""

# Authors: Chetan Gohil <chetan.gohil@psych.ox.ac.uk>

import os
import pathlib
from glob import glob
from dask.distributed import Client
from osl import source_recon, utils

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

# Directories
BASE_DIR = "/media/clement/CamCAN/MEG"
PREPROC_DIR = BASE_DIR + "/derivatives/preproc"
COREG_DIR = BASE_DIR + "/derivatives/coreg"
ANAT_DIR = "/media/clement/CamCAN/MRI/cc700-mri/mri/pipeline/release004/BIDS_20190411/anat"
FSL_DIR = "/home/clement/Bureau/fsl"
SRC_DIR = BASE_DIR + "/derivatives/src"

# Files
PREPROC_FILE = (
    PREPROC_DIR
    + "/mf2pt2_{0}_ses-rest_task-rest_meg"
    + "/mf2pt2_{0}_ses-rest_task-rest_meg_preproc_raw.fif"
)
SMRI_FILE = ANAT_DIR + "/{0}/anat/{0}_T1w.nii.gz"

config = """
    source_recon:
    - forward_model:
        model: Single Layer
    - beamform_and_parcellate:
        freq_range: [1, 45]
        chantypes: [mag, grad]
        rank: {meg: 60}
        parcellation_file: Glasser52_binary_space-MNI152NLin6_res-8x8x8.nii.gz
        method: spatial_basis
        orthogonalisation: symmetric
"""

if __name__ == "__main__":
    utils.logger.set_up(level="INFO")
    source_recon.setup_fsl(FSL_DIR)

    # Copy coreg directory
    # if not os.path.exists(SRC_DIR):
    #     cmd = f"cp -r {COREG_DIR} {SRC_DIR}"
    #     print(cmd)
    #     os.system(cmd)

    # Get subjects
    subjects = []
    for subject in sorted(
        glob(
            PREPROC_DIR
            + "/mf2pt2_*_ses-rest_task-rest_meg"
            + "/mf2pt2_sub-*_ses-rest_task-rest_meg_preproc_raw.fif"
        )
    ):
        subjects.append(pathlib.Path(subject).stem.split("_")[1])

    # Setup files
    smri_files = []
    preproc_files = []
    for subject in subjects:
        smri_files.append(SMRI_FILE.format(subject))
        preproc_files.append(PREPROC_FILE.format(subject))

    # Setup parallel processing
    client = Client(n_workers=8, threads_per_worker=1)

    # Run beamforming and parcellation
    source_recon.run_src_batch(
        config,
        src_dir=SRC_DIR,
        subjects=subjects,
        preproc_files=preproc_files,
        smri_files=smri_files,
        dask_client=False,
    )