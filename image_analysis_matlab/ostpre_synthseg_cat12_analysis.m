clear
close all
if ~isunix
    datadir = 'C:\Users\justoh\Data\ostpre_bids\derivatives';
else 
    datadir = '/research/work/justoh/ostpre_bids/bids/CAT12.8.1_roi';
end

% figdir = 'C:\Users\justoh\OneDrive - University of Eastern Finland\jussi\docs\ostpre-brain\figures';
figdir = 'C:\Users\justoh\OneDrive - University of Eastern Finland\jussi\docs\ostpre-brain\figures68';
xlsfn_ss = fullfile(datadir,'ostpre_synthseg_summary_measures.xlsx');
xlsfn_cat = fullfile(datadir,'ostpre_cat12_summary_measures.xlsx');


% must match the subjects first!

Tss = readtable(xlsfn_ss);
Tcat = readtable(xlsfn_cat);

subj_ss = cellfun(@(s) str2num(s(5:end)),Tss{:,1},"UniformOutput",true);
ses_ss = cellfun(@(s) str2num(s(5:end)),Tss{:,2},"UniformOutput",true);
subj_cat = cellfun(@(s) str2num(s(5:end)),Tcat{:,1},"UniformOutput",true);
ses_cat = cellfun(@(s) str2num(s(5:end)),Tcat{:,2},"UniformOutput",true);
[idx1,idx2] = find_same_indexes_nihpd([subj_ss ses_ss],[subj_cat ses_cat]);
% reorder tables 
Tss = Tss(idx1,:);
Tcat = Tcat(idx2,:);


% acceptable_rows = Tcat.quality_percent > 72; % this is the primary threshold
acceptable_rows = Tcat.quality_percent > 68;  % this is what Reijo has used 

regions = {'Jack_signature_CT','hippocampus','entorhinal_thickness','ventricle','GM_volume','TIV'};

for i = 1:length(regions)
    aaa{1}(:,i) = Tcat.(regions{i})(acceptable_rows);
    aaa{2}(:,i) = Tss.(regions{i})(acceptable_rows)/1000;
    mini = min([aaa{1}(:,i); aaa{2}(:,i)]);
    maxi = max([aaa{1}(:,i); aaa{2}(:,i)]);
    figure('Name',regions{i}) 
    scatter(aaa{1}(:,i),aaa{2}(:,i),"filled")
    c(i) = corr(aaa{1}(:,i),aaa{2}(:,i));
    txt = [regions{i},' Correlation ',num2str(c(i))];
    title(txt);
    xlabel('CAT12 Volume/mL')
    ylabel('Synthseg Volume/mL')
    hold on
    plot([mini maxi],[mini maxi],'k', 'LineWidth',3)
    print(fullfile(figdir,regions{i}),'-dpng')
    hold off
    [brob,stats] = robustfit(aaa{1}(:,i),aaa{2}(:,i));
    outliers_ind = abs(stats.resid)>stats.mad_s*8;
    inliers_ind = ~outliers_ind;
    Toutliers{i} = Tss(outliers_ind,1:2);
    figure("Name",strcat(regions{i},'_robust'));
    scatter(aaa{1}(inliers_ind,i),aaa{2}(inliers_ind,i),"blue","filled")
    cr(i) = corr(aaa{1}(inliers_ind,i),aaa{2}(inliers_ind,i));
    hold on
    bias(i) = mean(aaa{1}(inliers_ind,i) - aaa{2}(inliers_ind,i));
    scatter(aaa{1}(outliers_ind,i),aaa{2}(outliers_ind,i),"red","filled")
    txt = [regions{i}, ' RobCorr ',num2str(cr(i),2), ' Bias ', num2str(bias(i),4) , ' mL'];
    title(txt);
    xlabel('CAT12 Volume/mL')
    ylabel('Synthseg Volume/mL')
    plot([mini maxi],[mini maxi],'k', 'LineWidth',3)
    print(fullfile(figdir,strcat(regions{i},'_robust')),'-dpng')
    hold off
    
end




