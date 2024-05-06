
% Copyright 2023 Jussi Tohka, University of Eastern Finland jussi.tohka <at> uef.fi  
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
% documentation files (the �Software�), to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and 
% to permit persons to whom the Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED �AS IS�, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
% THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

clear
if ~isunix
    datadir = '\\research.uefad.uef.fi\groups\ostpre\ostpre_bids\bids\derivatives\CAT12.8.1';
    savedir = 'C:\Users\justoh\Data\ostpre_bids\derivatives';
else 
    datadir = '/research/groups/ostpre/ostpre_bids/bids/derivatives/CAT12.8.1';
    savedir = '/research/work/justoh/ostpre_bids/bids/CAT12.8.1_roi';
    atlasdir = '/research/users/justoh/matlab/spm12/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym';
end

atlasname ='aparc_DK40';
surface = 's'; % give surface = 's'; if surface
dd = dir(fullfile(datadir,'sub*'));
iter = 1;
for i = 1:length(dd)
    disp(i)
    ddd = dir(fullfile(datadir,dd(i).name,'ses*'));
    for j = 1:length(ddd)
       try 
           load(fullfile(datadir,dd(i).name,ddd(j).name,'anat',strcat('cat_',dd(i).name,'_',ddd(j).name,'_T1w.mat')));  % this is for TIV
           tab.vol_tiv(iter) = S.subjectmeasures.vol_TIV;  % TIV in ml
           tab.surf_TSA(iter) = S.subjectmeasures.surf_TSA;     % 
           tab.vol_abs_WMH(iter) = S.subjectmeasures.vol_abs_WMH; % white matter hyperintensities
           tab.IQR(iter)         = S.qualityratings.IQR; % quality ratings
           % turned to percentages by 
           % mark2rps    = @(mark) min(100,max(0,105 - real(mark)*10)) + isnan(real(mark)).*real(mark); 
           tab.tissue_typt_volCGW{iter} =  S.subjectmeasures.vol_abs_CGW; % these are in order CSF, gray matter, white matter according to cat12
           tab.average_cortical_thickness{iter} = S.subjectmeasures.dist_thickness{1}; % a 2 element vector, 1st component mean and 2nd std. The measure is very similar 
                                                                                 % to the one in the log files, but not exactly the same  
           load(fullfile(datadir,dd(i).name,ddd(j).name,'anat',strcat('catROI',surface,'_',dd(i).name,'_',ddd(j).name,'_T1w.mat')));
           tab.roithickness(iter,:) =  getfield(S,atlasname,'data','thickness')';
           tab.subj{iter} = dd(i).name;
           tab.session{iter} = ddd(j).name; 
           iter = iter + 1;
       catch
           disp(['Not found: ', fullfile(datadir,dd(i).name,ddd(j).name,'anat',strcat('cat_',dd(i).name,'_',ddd(j).name,'_T1w.mat'))]);
       end
    end
end
tab.roiids = getfield(S,atlasname,'ids');
tab.roinames = getfield(S,atlasname,'names');
% load(fullfile(datadir,dd(i).name,ddd(j).name,'anat',strcat('cat_',dd(i).name,'_',ddd(j).name,'_T1w.mat')));
% tab.software = S.software;
ipr = min(100,max(0,105 - real(tab.IQR)*10)) + isnan(real(tab.IQR)).*real(tab.IQR);

save(fullfile(savedir,strcat('roi_',atlasname)),'tab');



