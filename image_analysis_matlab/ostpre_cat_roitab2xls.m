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

clear
if ~isunix
    datadir = 'C:\Users\justoh\Data\ostpre_bids\derivatives';
else 
    datadir = '/research/work/justoh/ostpre_bids/bids/CAT12.8.1_roi';
end

atlas{1} ='roi_neuromorphometrics.mat';
atlas{2}  = 'roi_aparc_DK40.mat';
f(1) = load(fullfile(datadir,atlas{1}));
f(2) = load(fullfile(datadir,atlas{2}));

regions{1} = {'lentorhinal','rentorhinal','lfusiform','rfusiform','linferiortemporal','rinferiortemporal','lmiddletemporal','rmiddletemporal'};
regions{2} = {'Right Hippocampus','Left Hippocampus'};
regions{3} = {'lentorhinal','rentorhinal'};
regions{4} = {'Right Inf Lat Vent','Left Inf Lat Vent','Right Lateral Ventricle','Left Lateral Ventricle'};  

atlastype = [2 1 2 1];

nsubj = length(f(1).tab.subj);
varTypes = ["string","string","double","double","double","double","double","double","double","double"];
varNames = ["Subject","Session","Quality","TIV","mean CT","GM volume","Jack summary CT","Hippocampus volume","Entorhinal CT","Ventricle volume"];
sz = [nsubj,length(varTypes)];

S.subj = f(1).tab.subj;
S.session = f(1).tab.session;
S.vol_TIV = f(1).tab.vol_TIV;
T = struct2table(S);


if 0
ttt{1} = f(1).tab.roivoltot;
ttt{2} = f(2).tab.roithickness; 

for i = 1:length(regions)
    idx = zeros(length(regions{i}),1);
    for j = 1:length(regions{i})
        idx(j) = find(strcmp(regions{i}{j},f(atlastype(i)).tab.roinames));
    end
    sss{i} = ttt{atlastype(i)}(:,idx);
end
end

