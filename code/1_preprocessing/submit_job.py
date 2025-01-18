"""Submit jobs to the HPC cluster.

"""

import os
import numpy as np

def write_job_script(n_states, run):
    with open("job.sh", "w") as file:
        name = f"hmm_{n_states}_{run}"
        file.write("#!/bin/bash\n\n")
        file.write(f"#OAR -n {name}\n")
        file.write("#OAR -l /nodes=1/gpu=1,walltime=07:00:00\n")
        file.write(f"#OAR --stdout logs/{name}.out\n")
        file.write(f"#OAR --stderr logs/{name}.err\n")
        file.write("#OAR --project pr-deepneuro\n")
        file.write("#OAR -p gpumodel='A100'\n\n")

        file.write("cd /silenus/PROJECTS/pr-deepneuro/guichetc/code/\n")
        file.write("source /applis/environments/cuda_env.sh bigfoot 11.2\n")
        file.write("source /applis/environments/conda.sh\n")
        file.write("conda activate osld\n\n")
        file.write(f"python 01_TDE_HMM.py {n_states} {run}\n")

os.makedirs("logs", exist_ok=True)

for n_states in np.arange(6,12+2,2):
    for run in np.arange(1,5+1,1):
        write_job_script(n_states, run)
        name = f"hmm_{n_states}_{run}"
        os.system(f"mv job.sh ./{name}_job.sh")
        os.system(f"chmod +x ./{name}_job.sh")
        os.system(f"oarsub -S ./{name}_job.sh")

os.system("oarstat -u")