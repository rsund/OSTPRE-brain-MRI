% A script to traverse into CAT processed OSTPRE data and count the number of
% files in each sub-directory
clear
if isunix
    datadir = '/research/groups/ostpre/ostpre_bids/bids';
else
    datadir = '\\research.uefad.uef.fi\groups\ostpre/ostpre_bids/bids';
end
processed_dir_long = fullfile(datadir,'derivatives','CAT12.8.1_long');
processed_dir_cross = fullfile(datadir,'derivatives','CAT12.8.1');

subjid_fn = fullfile(datadir,'subjid.txt'); % note that this needs to be current
subjid_fid = fopen(subjid_fn,'r');
subjid = textscan(subjid_fid,'%s');
fclose(subjid_fid);

ns = length(subjid{1});

% d_long = dir(fullfile(processed_dir_long,'sub*'));
d_cross = dir(fullfile(processed_dir_cross,'sub*'));

% for each subject processed check the number of processing order 
% and then count the sessions and then files within the sessions
% store this info to the subject_status array 
% subject_status(:,1) = number_of_processed_sessions
% subject_status(:,2) = number_of_fully_processed_sessions
% fully processed
% should have 30 files, change 400623 by JT, reflected in reprocess_cross2.txt file   
% subject_status(:,3) = number_of_files_in_session_1
% subject_status(:,4) = number_of_files_in_session_2


% subject_status = zeros(ns,12);
% for i = 1:length(d_long)
%     s1 = d_long(i).name;
%     idx = find(strcmp(s1,subjid{1}));
%     dtmp = dir(fullfile(processed_dir_long,d_long(i).name,'ses*'));
%     subject_status(idx,1)= length(dtmp);
%     for j = 1:min(length(dtmp),11)
%         subject_status(idx,j + 1) = length(dir(fullfile(processed_dir_long,d_long(i).name,dtmp(j).name,'anat')));
%     end
% end

subject_status_c = zeros(ns,12);
for i = 1:length(d_cross)
    s1 = d_cross(i).name;
    idx = find(strcmp(s1,subjid{1}));
    dtmp = dir(fullfile(processed_dir_cross,d_cross(i).name,'ses*'));
    subject_status_c(idx,1)= length(dtmp);
    subject_status_c(idx,2) = 0;
    for j = 1:min(length(dtmp),11)
        subject_status_c(idx,j + 2) = length(dir(fullfile(processed_dir_cross,d_cross(i).name,dtmp(j).name,'anat'))); % change 040623 JT
        if length(dir(fullfile(processed_dir_cross,d_cross(i).name,dtmp(j).name,'anat'))) > 29
           subject_status_c(idx,2) = subject_status_c(idx,2) + 1;
        end
    end
    
end 
anat_files_per_subj = zeros(ns,1);
for i = 1:length(subjid{1})
    count = 0;
    dtmp = dir(fullfile(datadir,subjid{1}{i},'ses*'));
    for j = 1:length(dtmp)
        count = count + (length(dir(fullfile(datadir,subjid{1}{i},dtmp(j).name,'anat','*T1w*'))) > 0);
    end
    anat_files_per_subj(i) = count;
end
noanat = find(anat_files_per_subj == 0);
pp = find(subject_status_c(:,1) == 0);
to_be_processed = setdiff(pp,noanat);        
% images_missing = anat_files_per_subj - subject_status(:,1);
% images_missing = images_missing.*(subject_status(:,1) > 1); % if just one session, then long should not exist
images_missing_c = anat_files_per_subj - subject_status_c(:,2); % change 040623 by JT 
% save C:\Users\justoh\Results\ADNI_processing\adni_processing_status251022 
save C:\Users\justoh\Results\ostpre_processing\ostpre_processing_statsus070923

if 1
indc = find(images_missing_c);
% fid = fopen(fullfile(datadir,'reprocess_cross3.txt'),'w');
for i = 1:length(indc)
  %  fprintf(fid,'%s\n',subjid{1}{indc(i)});
    cccc{i} =  subjid{1}{indc(i)};
end
% fclose(fid);
end

% indc = find(images_missing);
% fid = fopen(fullfile(datadir,'reprocess_long.txt'),'w');
% for i = 1:length(indc)
%     fprintf(fid,'%s\n',subjid{1}{indc(i)});
% end
% fclose(fid);
