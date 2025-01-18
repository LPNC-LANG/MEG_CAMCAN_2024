clearvars
clc

%% CON
vol = spm_vol('E:\Research_Projects\MEG_CamCAN\network_masks\networks\CAB-NP_volumetric_liberal.nii');
[CON, ~] = spm_read_vols(vol);
CON(CON~=4)=0; % CON
CON = double(CON~=0); % binarize
vol_new = vol;
vol_new.fname = 'CONbin.nii';
spm_write_vol(vol_new, CON);

%% DMN
vol = spm_vol('E:\Research_Projects\MEG_CamCAN\network_masks\networks\CAB-NP_volumetric_liberal.nii');
[DMN, ~] = spm_read_vols(vol);
DMN(DMN~=9)=0; % DMN
DMN = double(DMN~=0); % binarize
vol_new = vol;
vol_new.fname = 'DMNbin.nii';
spm_write_vol(vol_new, DMN);

%% FPN
vol = spm_vol('E:\Research_Projects\MEG_CamCAN\network_masks\networks\CAB-NP_volumetric_liberal.nii');
[FPN, ~] = spm_read_vols(vol);
FPN(FPN~=7)=0; % FPN
FPN = double(FPN~=0); % binarize
vol_new = vol;
vol_new.fname = 'FPNbin.nii';
spm_write_vol(vol_new, FPN);


%% LANG
vol = spm_vol('E:\Research_Projects\MEG_CamCAN\network_masks\networks\CAB-NP_volumetric_liberal.nii');
[LANG, ~] = spm_read_vols(vol);
LANG(LANG~=6)=0; % LANG
LANG = double(LANG~=0); % binarize
vol_new = vol;
vol_new.fname = 'LANGbin.nii';
spm_write_vol(vol_new, LANG);

%% SMN
vol = spm_vol('E:\Research_Projects\MEG_CamCAN\network_masks\networks\CAB-NP_volumetric_liberal.nii');
[SMN, ~] = spm_read_vols(vol);
SMN(SMN~=3)=0; % SMN
SMN = double(SMN~=0); % binarize
vol_new = vol;
vol_new.fname = 'SMNbin.nii';
spm_write_vol(vol_new, SMN);
%% SCN
vol = spm_vol('E:\Research_Projects\MEG_CamCAN\network_masks\networks\SCN.nii');
[SCN, ~] = spm_read_vols(vol);
SCN = double(SCN~=0); % binarize
vol_new = vol;
vol_new.dt = [spm_type('float32') 0];
vol_new.fname = 'SCNbin.nii';
spm_write_vol(vol_new, SCN);

%% MD
vol = spm_vol('E:\Research_Projects\MEG_CamCAN\network_masks\networks\MD_Fedorenko.nii.gz');
[MD, ~] = spm_read_vols(vol);
MD = double(MD~=0); % binarize
vol_new = vol;
vol_new.fname = 'MDbin.nii';
spm_write_vol(vol_new, MD);