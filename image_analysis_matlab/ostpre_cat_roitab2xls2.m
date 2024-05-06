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
xlsfn = fullfile(datadir,'ostpre_CAT12_summary_measures.xlsx');
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

ct = cell2mat(f(1).tab.average_cortical_thickness');
totvol = cell2mat(f(1).tab.tissue_typt_volCGW');

S.subj = f(1).tab.subj';
S.session = f(1).tab.session';
S.quality_percent = (min(100,max(0,105 - real(f(1).tab.IQR)*10)) + isnan(real(f(1).tab.IQR)).*real(f(1).tab.IQR))';
S.quality_IQR = f(1).tab.IQR';
S.TIV = f(1).tab.vol_tiv';
S.average_thickness = ct(:,1);
S.GM_volume = totvol(:,2);


ttt{1} = f(1).tab.roivoltot;
ttt{2} = f(2).tab.roithickness; 

for i = 1:length(regions)
    idx = zeros(length(regions{i}),1);
    for j = 1:length(regions{i})
        idx(j) = find(strcmp(regions{i}{j},f(atlastype(i)).tab.roinames));
    end
    sss{i} = ttt{atlastype(i)}(:,idx);

end
S.Jack_signature_CT = mean(sss{1},2);
S.hippocampus = sum(sss{2},2);
S.entorhinal_thickness = mean(sss{3},2);
S.ventricle = sum(sss{4},2);
T = struct2table(S);
writetable(T,xlsfn);
