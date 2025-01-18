% Clément Guichet, UGA CNRS UMR 5105 LPNC, May 2024
%% IMPORT DATA

clc
clearvars
close all

n_states = 8;
% Temporal data from the state time courses of each subject
fo = readNPY("E:/Research_Projects/MEG_CamCAN/TDE_HMM/results/inf_params/08_states/fo_reordered.npy");
lt = readNPY("E:/Research_Projects/MEG_CamCAN/TDE_HMM/results/inf_params/08_states/lt_reordered.npy");
intv = readNPY("E:/Research_Projects/MEG_CamCAN/TDE_HMM/results/inf_params/08_states/intv_reordered.npy");
sr = readNPY("E:/Research_Projects/MEG_CamCAN/TDE_HMM/results/inf_params/08_states/sr_reordered.npy");

% Concatenating and zscoring across states
brain_matrix = cat(2, ...
    fo,lt,intv,sr ...
    );

disp("Discarded subjects with missing >3 cognitive variables");
disp("Subjects selected");
mask = importdata("E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis/mask_subjects_education.csv").data;
brain_matrix = brain_matrix(logical(mask),:); % Keep subjects with cog data

brain_labels = cell(4*n_states,1);
for state = 1:n_states; brain_labels{state,1} = ['fo ' num2str(state)]; end;
for state = 1:n_states; brain_labels{n_states + state,1} = ['lt ' num2str(state)]; end;
for state = 1:n_states; brain_labels{2*n_states + state,1} = ['intv ' num2str(state)]; end;
for state = 1:n_states; brain_labels{3*n_states + state,1} = ['sr ' num2str(state)]; end;

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
disp("*Quantile normalization - DONE*");
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
quadratic_term = -shifted_age.^2;

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
%     quadratic_cog_age(idx) = quadratic_term(value-18);
    level_off_cog_age(idx) = level_off_term(value-18);
    acceleration_cog_age(idx) = acceleration_term(value-18);
end

figure
hold on
% scatter(cog_age,cog_age)
% plot(cog_age,quadratic_cog_age)
plot(cog_age,level_off_cog_age)
plot(cog_age,acceleration_cog_age)
hold off
%% Heatmap for Figure 1
H = heatmap(figure,cog_matrix_preprocessed,...
    'ColorScaling','scaledrows',...
    'FontName','Arial',...
    'CellLabelFormat','%0.2g',...
    'XDisplayLabels',cog_data.Properties.VariableNames(1,7:15));

hHeatmap = struct(H).Heatmap;
hHeatmap.GridLineStyle = ':';

%%
cog_matrix_age = cat(2,level_off_cog_age.', acceleration_cog_age.', quantilenorm(cog_matrix_preprocessed));
cog_labels = cat(2,...
    "level_off", "accelerate", ...
    cog_data.Properties.VariableNames(1,7:15)...
    );

disp("*Data preparation - DONE*");
clearvars -except cog_matrix_age cog_labels ...
    brain_matrix brain_labels ...
    grouping n_states

%% Check all inputs for validity
grouped_plots = 0;

output_path = 'E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis/output_PLS_temporal'; 
inputs_PLS_temporal
[input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);
%% Run PLS analysis (including permutation testing and bootstrapping)
res = myPLS_analysis(input,pls_opts); 

%% BSR values

BSR_U = res.U./res.boot_results.Ub_std;
BSR_age = BSR_U(1:2,:);
BSR_cog = BSR_U(3:end,:);
BSR_ratios_brain = res.V./res.boot_results.Vb_std;
%% Save results struct
disp('... Saving results ...');
save(fullfile(save_opts.output_path,[save_opts.prefix '_res']),'res','-v7.3');
disp(' ');
myPLS_plot_results(res,save_opts);
%%
% grouped_plots = 1;
% output_path = './output_PLS_temporal'; 
% inputs_PLS_temporal
% [input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);
% 
% 
% groupIDs=unique(grouping);
% nGroups=length(groupIDs);
% names_groups= {'18-44',...
%                    '45-55',...
%                    '56+'};
%                
% colors = {'b','r','c','g','m','y','w','k'}; % Matlab colors
% plot_colors = colors(1:nGroups); % select as many colors as groups
% 
% disp('Correlations between imaging and behavioral scores');
% for iter_lc = 1:2
%     this_lc = iter_lc;
%     
%     figure('position',[440   541   327   257]);
%     for iG = 1:nGroups
%         plot(res.Lx(grouping==groupIDs(iG),this_lc),...
%             res.Ly(grouping==groupIDs(iG),this_lc),[plot_colors{iG} '.'],'MarkerSize',20);
%         hold on
%     end
%     hold off
%     title(['LC' num2str(this_lc) ' - Correlations between imaging and behavioral scores']);
%     
%     
%     legend(names_groups,'Location','southeast');
%     xlabel('Imaging scores');
%     ylabel('Behavior/Design scores');
%     
%     set(gcf,'Color','w');
%     set(gca,'Box','off');
%     
%     corr_LxLy(iter_lc) = corr(res.Lx(:,this_lc),res.Ly(:,this_lc));
%     disp(['LC' num2str(this_lc) ': r = ' num2str(corr_LxLy(iter_lc),'%0.2f')]);
%     
%     print(gcf,'-depsc2','-painters');
% end
