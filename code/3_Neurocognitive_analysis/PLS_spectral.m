% Clément Guichet, UGA CNRS UMR 5105 LPNC, May 2024
%% IMPORT DATA

clc
clearvars
close all

%%
% Spectral data from the state time courses of each subject
% n_subjects = 597; 
n_subjects = 467; % when adding in education levels
n_states = 8;
n_parcels = 20;
n_freqs = 88;

brain_data = readNPY("E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis/psd_pls.npy");
size(brain_data);
%%
disp("Discarded subjects with missing >3 cognitive variables");
disp("Subjects selected");
mask = importdata("E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis/mask_subjects_education.csv").data;
brain_matrix = brain_data(logical(mask),:); % Keep subjects with cog data
brain_matrix_reshaped = reshape(brain_matrix, [n_subjects,n_states*n_parcels*n_freqs]);

brain_mask = logical(brain_matrix_reshaped); % Keep logical for indexing after PLS
disp("Logical mask generated");
brain_matrix_relevant = brain_matrix_reshaped(:,all(brain_matrix_reshaped~=0)); % Discard all zero entries
disp("Discarded frequency bins not relevant for each state");

brain_labels = cell(n_states*n_parcels*n_freqs,1);
brain_labels_relevant = brain_labels(all(brain_matrix~=0));
disp("Assigning brain labels");
%% Define all the inputs
% Cognitive data
cog_data = readtable("E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis/cog_data_education.csv");
% Cattell, Proverb, Naming, ToT_Ratio, Hotel_Task, Sentence_Comprehension, 
% Story_Recall, Verbal_Fluency, YA_edcuation 
cog_matrix = table2array(cog_data(:,7:15)); 

% Gender, TIV, MMSE, nsamples
nsamples = readNPY("E:/Research_Projects/MEG_CamCAN/TDE_HMM/results/inf_params/08_states/w.npy");
nsamples_final = nsamples(logical(mask),:);
covariates = cat(2, table2array(cog_data(:,4:6)), nsamples_final); 
% Regress out covariates
[cog_matrix_preprocessed, ~, ~, ~] = CBIG_glm_regress_matrix(...
    cog_matrix,...
    cat(2, covariates(:,1),zscore(covariates(:,2:end)))... % Gender contrast need not be normalized
    );
disp("*Regressing covariates out - DONE*");

%% Adding age into the matrix
cog_age = table2array(cog_data(:,3)); % Age

% Create piece-wise contrasts
age = (19:89)';
midlife_age = 55;

% Compute the linear term
linear_term = age;
% Compute the quadratic term centered at midlife
shifted_age = age - midlife_age;
quadratic_term = -shifted_age.^2; % Inverted U-shape

% Normalize the quadratic term to match the scale of the linear term
quadratic_term = (quadratic_term - min(quadratic_term)) / (max(quadratic_term) - min(quadratic_term));
quadratic_term = quadratic_term * (max(linear_term) - min(linear_term)) + min(linear_term);

% Compute the term that levels off at midlife
level_off_term = zeros(size(age));
level_off_term(age <= midlife_age) = quadratic_term(age <= midlife_age);
level_off_term(age > midlife_age) = quadratic_term(age == midlife_age);

% Mirror the level off behavior for the acceleration term
acceleration_term = zeros(size(age));
acceleration_term(age <= midlife_age) = quadratic_term(age == midlife_age);
acceleration_term(age > midlife_age) = quadratic_term(age > midlife_age);

%% Now apply the piecewise functions to the continuous age values in the sample
for idx = 1:size(cog_age) 
    value = round(cog_age(idx)); % Grab the age value to retrieve
    quadratic_cog_age(idx) = quadratic_term(value-18);
    level_off_cog_age(idx) = level_off_term(value-18);
    acceleration_cog_age(idx) = acceleration_term(value-18);
end

figure
hold on
scatter(cog_age,cog_age)
plot(cog_age,quadratic_cog_age)
plot(cog_age,level_off_cog_age)
plot(cog_age,acceleration_cog_age)
hold off
%%
cog_matrix_age = cat(2,level_off_cog_age.',acceleration_cog_age.',quantilenorm(cog_matrix_preprocessed));
cog_labels = cat(2,...
    "level_off", "accelerate", ...
    cog_data.Properties.VariableNames(1,7:15)...
    );

disp("*Quantile normalization - DONE*");
disp("*Data preparation - DONE*");
clearvars -except cog_matrix_age cog_labels ...
    brain_matrix_relevant brain_labels_relevant ...
    brain_mask ...
    n_states n_parcels n_freqs
%% Check all inputs for validity
grouped_plots = 0;

output_path = 'E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis/output_PLS_spectral'; 
inputs_PLS_spectral
[input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);
%% Run PLS analysis (including permutation testing and bootstrapping)
res = myPLS_analysis(input,pls_opts); 

%% Bootstrap sampling ratios
BSR_U = res.U./res.boot_results.Ub_std; % Cognition
BSR_age = BSR_U(1:2,:);
BSR_cog = BSR_U(3:end,:);

BSR_ratios_brain = res.V./res.boot_results.Vb_std; % Brain
LC1_BSR_brain = BSR_ratios_brain(:,1);

n_components = find(res.LC_pvals<= 0.05, 1, 'last'); % Find number of significant components

brain_mask = logical(brain_mask(1,:)); % Keep only one row for iterating

% Pre-allocation
BSR_expanded = zeros([1,size(brain_mask,2)]);
BSR_PLS_bigmodel = zeros([n_components, n_states, n_parcels, n_freqs]); 

for comp = 1:n_components
    BSR_tmp = BSR_ratios_brain(:,comp).'; % Grab BSR values of the component
    BSR_expanded(brain_mask) = BSR_tmp; % Constrain values to logical mask
    BSR_PLS_bigmodel(comp,:,:,:) = reshape(BSR_expanded, [1, n_states, n_parcels, n_freqs]);
end

writeNPY(BSR_PLS_bigmodel, ...
    fullfile(output_path,['/BSR_PLS_spectral.npy']));

%% Save results struct
disp('... Saving results ...');
save(fullfile(save_opts.output_path,[save_opts.prefix '_res']),'res','-v7.3');
disp(' ');
myPLS_plot_results(res,save_opts);

%%
myPLS_plot_subjScores(res.Lx,res.Ly,res.group_names,res.grouping,2,save_opts);