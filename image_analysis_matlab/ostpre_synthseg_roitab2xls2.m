% this produces a limited number of T1w MRI features based on 
% volumetric ROIs and cortical thickness values
% hippocampus
% ventricles
% entorhinal cortex
% total gray matter
% icv
% total thickness
% Jack et al summary measure
% A measure of cortical thickness 
% in AD signature regions was calculated as the average of cortical thickness in entorhinal,
%  inferior temporal, middle temporal, and fusiform regions as previously described 
% Jack CR Jr, Wiste HJ, Weigand SD, Knopman DS, Mielke MM, Vemuri P, Lowe V, Senjem ML, Gunter JL, Reyes D, 
% Machulda MM, Roberts R, Petersen RC (2015) Different definitions of neurodegeneration produce similar 
% amyloid/neurodegeneration biomarker group findings. Brain 138, 3747–3759. [
% This script is for synthseg 

clear
if ~isunix
    datadir = 'C:\Users\justoh\Data\ostpre_bids\derivatives';
else 
    datadir = '/research/work/justoh/ostpre_bids/bids/CAT12.8.1_roi';
end
xlsfn = fullfile(datadir,'ostpre_synthseg_summary_measures.xlsx');
atlas{1} ='roi_neuromorphometrics.mat';
atlas{2}  = 'roi_aparc_DK40.mat';
Tf = readtable('C:\Users\justoh\OneDrive - University of Eastern Finland\jussi\docs\ostpre-brain\synthseg_segmentations\synthseg_vol.xlsx');
% f(2) = load(fullfile(datadir,atlas{2}));

% note that variable names need to be changed for Matlab. readtable does
% that automatically
regions{1} = {'ctx_lh_entorhinal','ctx_rh_entorhinal','ctx_lh_fusiform','ctx_rh_fusiform','ctx_lh_inferiortemporal','ctx_rh_inferiortemporal','ctx_lh_middletemporal','ctx_rh_middletemporal'};
regions{2} = {'right_hippocampus','left_hippocampus'};
regions{3} = {'ctx_lh_entorhinal','ctx_rh_entorhinal'}; % note that these are volumes so not totally comparable to the cortical thickness ones
regions{4} = {'right_inferior_lateral_ventricle','left_inferior_lateral_ventricle','right_lateral_ventricle','left_lateral_ventricle'};  

% atlastype = [2 1 2 1];

nsubj = size(Tf,1);
varTypes = ["string","string","double","double","double","double","double","double","double","double"];
varNames = ["Subject","Session","Quality","TIV","mean CT","GM volume","Jack summary CT","Hippocampus volume","Entorhinal CT","Ventricle volume"];
sz = [nsubj,length(varTypes)];

% ct = cell2mat(f(1).tab.average_cortical_thickness');
% totvol = cell2mat(f(1).tab.tissue_typt_volCGW');

for i = 1:nsubj
    S.subj{i,1} = Tf.X{i}(1:9);
    S.ses{i,1} = Tf.X{i}(11:15);
    S.quality_percent(i,1) = -1; % not available at the moment 
    S.quality_IQR(i,1) = -1; % not available at the moment
    S.TIV(i,1) = Tf.total_intracranial(i); 
    S.average_thickness(i,1) = -1; % not available in synthseg
    S.GM_volume(i,1) = Tf.left_cerebral_cortex(i) + Tf.right_cerebral_cortex(i) + Tf.left_cerebellum_cortex(i) + Tf.right_cerebellum_cortex(i); 
end

for i = 1:length(regions)
    idx = zeros(length(regions{i}),1);
    for j = 1:length(regions{i})
        idx(j) = find(strcmp(regions{i}{j},Tf.Properties.VariableNames));
    end
    sss{i} = Tf{:,idx};

end
S.Jack_signature_CT = mean(sss{1},2);
S.hippocampus = sum(sss{2},2);
S.entorhinal_thickness = mean(sss{3},2);
S.ventricle = sum(sss{4},2);
T = struct2table(S);
writetable(T,xlsfn);
