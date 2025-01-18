# Regrid the networks to match that of the MEG atlas
mrgrid SCNbin.nii regrid -template /mnt/d/Analyses/Atlases/HCP/Glasser52_binary_space-MNI152Lin6_res-8x8x8.nii r_SCNbin.nii
mrgrid MDbin.nii regrid -template /mnt/d/Analyses/Atlases/HCP/Glasser52_binary_space-MNI152Lin6_res-8x8x8.nii r_MDbin.nii
mrgrid DMNbin.nii regrid -template /mnt/d/Analyses/Atlases/HCP/Glasser52_binary_space-MNI152Lin6_res-8x8x8.nii r_DMNbin.nii
mrgrid CONbin.nii regrid -template /mnt/d/Analyses/Atlases/HCP/Glasser52_binary_space-MNI152Lin6_res-8x8x8.nii r_CONbin.nii
mrgrid SMNbin.nii regrid -template /mnt/d/Analyses/Atlases/HCP/Glasser52_binary_space-MNI152Lin6_res-8x8x8.nii r_SMNbin.nii
mrgrid LANGbin.nii regrid -template /mnt/d/Analyses/Atlases/HCP/Glasser52_binary_space-MNI152Lin6_res-8x8x8.nii r_LANGbin.nii
mrgrid FPNbin.nii regrid -template /mnt/d/Analyses/Atlases/HCP/Glasser52_binary_space-MNI152Lin6_res-8x8x8.nii r_FPNbin.nii

# Smoothing with a Gaussian kernel with FWHM of 10 mm
mrfilter r_SCNbin.nii smooth -fwhm 10 sr_SCNbin.nii
mrfilter r_MDbin.nii smooth -fwhm 10 sr_MDbin.nii
mrfilter r_DMNbin.nii smooth -fwhm 10 sr_DMNbin.nii
mrfilter r_CONbin.nii smooth -fwhm 10 sr_CONbin.nii
mrfilter r_SMNbin.nii smooth -fwhm 10 sr_SMNbin.nii
mrfilter r_LANGbin.nii smooth -fwhm 10 sr_LANGbin.nii
mrfilter r_FPNbin.nii smooth -fwhm 10 sr_FPNbin.nii

cd ../glasser52
# ls . | grep .nii > region_list52.txt
# sed 's/.nii//g' region_list52.txt > region_list52.txt

for region in `cat region_list52.txt`; do
	mrstats ${region}.nii -ignorezero | grep [0] >> ../output/MD_overlap.tsv
	mrstats ${region}.nii -mask ../networks/r_MDbin.nii -ignorezero | grep [0] >> ../output/MD_overlap.tsv

	mrstats ${region}.nii -ignorezero | grep [0] >> ../output/SCN_overlap.tsv
	mrstats ${region}.nii -mask ../networks/r_SCNbin.nii -ignorezero | grep [0] >> ../output/SCN_overlap.tsv

	mrstats ${region}.nii -ignorezero | grep [0] >> ../output/DMN_overlap.tsv
	mrstats ${region}.nii -mask ../networks/r_DMNbin.nii -ignorezero | grep [0] >> ../output/DMN_overlap.tsv

	mrstats ${region}.nii -ignorezero | grep [0] >> ../output/CON_overlap.tsv
	mrstats ${region}.nii -mask ../networks/r_CONbin.nii -ignorezero | grep [0] >> ../output/CON_overlap.tsv

	mrstats ${region}.nii -ignorezero | grep [0] >> ../output/LANG_overlap.tsv
	mrstats ${region}.nii -mask ../networks/r_LANGbin.nii -ignorezero | grep [0] >> ../output/LANG_overlap.tsv

	mrstats ${region}.nii -ignorezero | grep [0] >> ../output/SMN_overlap.tsv
	mrstats ${region}.nii -mask ../networks/r_SMNbin.nii -ignorezero | grep [0] >> ../output/SMN_overlap.tsv

	mrstats ${region}.nii -ignorezero | grep [0] >> ../output/FPN_overlap.tsv
	mrstats ${region}.nii -mask ../networks/sr_FPNbin.nii -ignorezero | grep [0] >> ../output/FPN_overlap.tsv
done
