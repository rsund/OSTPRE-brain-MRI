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
% amyloid/neurodegeneration biomarker group findings. Brain 138, 3747–3759.

% these are then provided for listed ADNI subjects given in a matchfn file

clear
% this file gives the subject RIDs, presumably baseline
% Here are the RIDs for 6 samples with LMCI at DX_bl and Dementia diagnosis in DX at VISCODE of bl: 332, 995, 1154, 78,190,1226.
matchfn = 'C:\Users\justoh\OneDrive - University of Eastern Finland\jussi\docs\ostpre-brain\adni_match_sample\ADNI_whiteWomen.xlsx';
Trid = readtable(matchfn);
rid = Trid.RID;
idx_remove = find(ismember(rid,[332 995 1154 78 190 1226])); % 
datadir = 'C:\Users\justoh\Data\ADNI_bids\derivatives';
xlsfn = fullfile(datadir,'ADNI_CAT12_ostprematched_summary_measures_v2.xlsx');

% OLD CODE: PROVIDED ONLY ABOUT 65 % OF THE SCANNER mf'S PROBABLY DUE TO THE REQUIREMENT OF T1,T2, AND PD SCANS  
% scannerfn = fullfile(datadir,'ADNI_UCD_WMH_05_02_22_07Mar2024.csv'); % this file contains the information about the scanner manufacturers
% Tscanner = readtable(scannerfn);
% idx_bl = find(strcmp(Tscanner.VISCODE2,'scmri') | strcmp(Tscanner.VISCODE2,'sc') | strcmp(Tscanner.VISCODE2,'bl'));
% scanner_rid = Tscanner.RID(idx_bl);
% [idx_scanner1 idx_scanner2] = find_same_indexes_nihpd([scanner_rid zeros(size(scanner_rid))],[rid zeros(size(rid))]);

% NEW IDEA: Search from LONI image archieve by Scanner MF and download the
% image meta information as csv-file, there's a screenshot on how to do this  
% remove the participants that cannot be decided based on this

scannerfn{1} = fullfile(datadir,'siemens_idaSearch_3_07_2024.csv'); % this file contains the information about the scanner manufacturers
scannerfn{2} = fullfile(datadir,'philips_idaSearch_3_07_2024.csv');
scannerfn{3} = fullfile(datadir,'GE_idaSearch_3_07_2024.csv');

for i = 1:3
    Tscanner = readtable(scannerfn{i});
    ssubj = Tscanner.SubjectID;
    for j = 1:length(ssubj)
        scanner_rid{i}(j) = str2num(ssubj{j}(7:end));
    end
end
scanner = zeros(length(rid),3);
for j = 1:3
    scanner(:,j) = ismember(rid,scanner_rid{j});
end
% these are subjects whose scanner information cannot be deduced this way
% 
idx_uncertain = find(sum(scanner,2) ~= 1);
idx_remove = union(idx_remove,idx_uncertain);
[~, scanner1d] = max(scanner,[],2);
idx_keep = setdiff(1:length(rid),idx_remove);
rid = rid(idx_keep);

atlas{1} ='roi_neuromorphometrics_v2.mat';
atlas{2}  = 'roi_aparc_DK40.mat';
f(1) = load(fullfile(datadir,atlas{1}));
f(2) = load(fullfile(datadir,atlas{2}));

regions{1} = {'lentorhinal','rentorhinal','lfusiform','rfusiform','linferiortemporal','rinferiortemporal','lmiddletemporal','rmiddletemporal'};
regions{2} = {'Right Hippocampus','Left Hippocampus'};
regions{3} = {'lentorhinal','rentorhinal'};
regions{4} = {'Right Inf Lat Vent','Left Inf Lat Vent','Right Lateral Ventricle','Left Lateral Ventricle'};  

atlastype = [2 1 2 1];

nsubj = length(f(1).tab.subj);

ct = cell2mat(f(1).tab.average_cortical_thickness');
totvol = cell2mat(f(1).tab.tissue_typt_volCGW');

for i = 1:nsubj
    fsubj(i) = str2num(f(1).tab.subj{i}(13:end));
    fses(i) = str2num(f(1).tab.session{i}(6:end));
end
[idx1 idx2] = find_same_indexes_nihpd([fsubj' fses'],[rid zeros(size(rid))]);
% note that idx2 referes to the indexes after uncertain participants have
% been taken away 

S.subj = rid(idx2);
S.session = f(1).tab.session(idx1)';
tmp = Trid.DX_bl(idx_keep);
S.DX_bl = tmp(idx2);
tmp = Trid.AGE(idx_keep);
S.Age = tmp(idx2);
tmp = Trid.FLDSTRENG(idx_keep);
S.FLDSTRENG = tmp(idx2);
tmp = Trid.ORIGPROT(idx_keep);
S.ORIGPROT = tmp(idx2);
for i = 1:length(S.FLDSTRENG)
   if isempty(S.FLDSTRENG{i})
       if strcmp(S.ORIGPROT{i},'ADNI1')
           S.FLDSTRENG{i} = '1.5 Tesla MRI';
       else
           S.FLDSTRENG{i} =  '3 Tesla MRI';
       end
   end
end
tmp = scanner1d(idx_keep);
S.MANUFACTURER = tmp(idx2);

S.quality_percent = (min(100,max(0,105 - real(f(1).tab.IQR)*10)) + isnan(real(f(1).tab.IQR)).*real(f(1).tab.IQR))';
S.quality_percent = S.quality_percent(idx1);
S.quality_IQR = f(1).tab.IQR(idx1)';
S.TIV = f(1).tab.vol_tiv(idx1)';
S.average_thickness = ct(idx1,1);
S.GM_volume = totvol(idx1,2);


ttt{1} = f(1).tab.roivoltot;
ttt{2} = f(2).tab.roithickness; 

for i = 1:length(regions)
    idx = zeros(length(regions{i}),1);
    for j = 1:length(regions{i})
        idx(j) = find(strcmp(regions{i}{j},f(atlastype(i)).tab.roinames));
    end
    sss{i} = ttt{atlastype(i)}(idx1,idx);

end
S.Jack_signature_CT = mean(sss{1},2);
S.hippocampus = sum(sss{2},2);
S.entorhinal_thickness = mean(sss{3},2);
S.ventricle = sum(sss{4},2);
T = struct2table(S);
writetable(T,xlsfn);