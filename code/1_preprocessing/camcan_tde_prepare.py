"""Data preparation: time-delay embedding and principal component analysis.

"""

# Authors: Chetan Gohil <chetan.gohil@psych.ox.ac.uk>

from glob import glob

from osl_dynamics.data import Data

files = sorted(glob("/media/clement/CamCAN/MEG/derivatives/src/*/sflip_parc-raw.fif"))
data = Data(files, picks="misc", reject_by_annotation="omit", n_jobs=8)
methods = {
    "tde_pca": {"n_embeddings": 15, "n_pca_components": 104},
    "standardize": {},
}
data.prepare(methods)
data.save("/media/clement/CamCAN/MEG/derivatives/HMM_prepared")
data.delete_dir()
