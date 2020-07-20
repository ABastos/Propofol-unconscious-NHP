

cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')


filestoload = ls('*_TFR_spikes*')

time_axis_global = -80:10/60:180;
TFR_prePropofol_bs_sessions = NaN(size(filestoload,1), 4,50,length(-20*60:10:-15*60));
TFRbs_sessions = NaN(size(filestoload,1), 4,50,length(-30*60:10:9*60));
TFR_sessions = NaN(size(filestoload,1), 4,50,length(-30*60:10:9*60));
drugstart_sessions = zeros(size(filestoload,1),1);
drugstart2_sessions = zeros(size(filestoload,1),1);
drugstop_sessions = zeros(size(filestoload,1),1);
LOC_sessions = zeros(size(filestoload,1),1);
ROC_sessions = zeros(size(filestoload,1),1);
for f = 1:size(filestoload,1)

    load(filestoload(f,:)); %, 'sessionInfo', 'spike_rates', 'time_axis_baseline')
    
    if isempty(sessionInfo.eyesOpen)
        error('deal with no eye info')
    end
    
    if isempty(sessionInfo.eyesClose)
        error('deal with no eye info')
    else
        time_axis_LOC = time_axis_baseline-((sessionInfo.eyesClose(1)-sessionInfo.drugStart(1))/60);
    end
   
    
    tbin1PP = nearest(time_axis_LOC, -20);
    tbin2PP = nearest(time_axis_LOC, -15);

    if f<=11
        time_axis_ROC = TFR.time-sessionInfo.eyesOpen(1);
        TFRbinstart = nearest(TFR.time, sessionInfo.drugStart(1));
        TFRbinstart2 = nearest(TFR.time, sessionInfo.drugStart(2));
        TFRbinstop = nearest(TFR.time, sessionInfo.drugEnd(2));
        TFRbinLOC = nearest(TFR.time, sessionInfo.eyesClose(1));
        TFRbinROC = nearest(TFR.time, sessionInfo.eyesOpen(1));
        
    tbin1 = nearest(time_axis_ROC, -30*60);
    tbin2 = nearest(time_axis_ROC, 9*60);        
        
        PFC_indx = [];
        FEF_indx = [];
        PPC_indx = [];
        STG_indx = [];
        for c = 1:length(TFRbs.label)
            if strcmp(TFRbs.label{c}(1:2), '7b')
                PPC_indx = [PPC_indx; c];
            elseif strcmp(TFRbs.label{c}(1:3), 'FEF')
                FEF_indx = [FEF_indx; c];
            elseif strcmp(TFRbs.label{c}(1:3), 'CPB')
                STG_indx = [STG_indx; c];
            elseif strcmp(TFRbs.label{c}(1:5), 'vlPFC')
                PFC_indx = [PFC_indx; c];
            end
        end
        
        TFRbs_sessions(f,1,:,:) = squeeze(squeeze(mean(TFRbs.powspctrm(1,PFC_indx,:,tbin1:tbin2),2)));
        TFRbs_sessions(f,2,:,:) = squeeze(squeeze(mean(TFRbs.powspctrm(1,FEF_indx,:,tbin1:tbin2),2)));
        TFRbs_sessions(f,3,:,:) = squeeze(squeeze(mean(TFRbs.powspctrm(1,PPC_indx,:,tbin1:tbin2),2)));
        TFRbs_sessions(f,4,:,:) = squeeze(squeeze(mean(TFRbs.powspctrm(1,STG_indx,:,tbin1:tbin2),2)));
        %Absolute
        TFR_sessions(f,1,:,:) = squeeze(squeeze(mean(TFR.powspctrm(1,PFC_indx,:,tbin1:tbin2),2)));
        TFR_sessions(f,2,:,:) = squeeze(squeeze(mean(TFR.powspctrm(1,FEF_indx,:,tbin1:tbin2),2)));
        TFR_sessions(f,3,:,:) = squeeze(squeeze(mean(TFR.powspctrm(1,PPC_indx,:,tbin1:tbin2),2)));
        TFR_sessions(f,4,:,:) = squeeze(squeeze(mean(TFR.powspctrm(1,STG_indx,:,tbin1:tbin2),2))); 
        %Pre-propofol
        TFR_prePropofol_bs_sessions(f,1,:,:) = squeeze(squeeze(mean(TFRbs.powspctrm(1,PFC_indx,:,tbin1PP:tbin2PP),2)));
        TFR_prePropofol_bs_sessions(f,2,:,:) = squeeze(squeeze(mean(TFRbs.powspctrm(1,FEF_indx,:,tbin1PP:tbin2PP),2)));
        TFR_prePropofol_bs_sessions(f,3,:,:) = squeeze(squeeze(mean(TFRbs.powspctrm(1,PPC_indx,:,tbin1PP:tbin2PP),2)));
        TFR_prePropofol_bs_sessions(f,4,:,:) = squeeze(squeeze(mean(TFRbs.powspctrm(1,STG_indx,:,tbin1PP:tbin2PP),2)));
    else
        time_axis_ROC = ALL_TFR{1}.time-sessionInfo.eyesOpen(1);
        
        TFRbinstart = nearest(ALL_TFR{1}.time, sessionInfo.drugStart(1));
        TFRbinstart2 = nearest(ALL_TFR{1}.time, sessionInfo.drugStart(2));
        TFRbinstop = nearest(ALL_TFR{1}.time, sessionInfo.drugEnd(2));
        TFRbinLOC = nearest(ALL_TFR{1}.time, sessionInfo.eyesClose(1));
        TFRbinROC = nearest(ALL_TFR{1}.time, sessionInfo.eyesOpen(1));
        
        
        TFRbs_sessions(f,1,:,:) = squeeze(squeeze(mean(ALL_TFRbs{1}.powspctrm(1,:,:,tbin1:tbin2),2)));
        TFRbs_sessions(f,2,:,:) = squeeze(squeeze(mean(ALL_TFRbs{2}.powspctrm(1,:,:,tbin1:tbin2),2)));
        TFRbs_sessions(f,3,:,:) = squeeze(squeeze(mean(ALL_TFRbs{3}.powspctrm(1,:,:,tbin1:tbin2),2)));
        TFRbs_sessions(f,4,:,:) = squeeze(squeeze(mean(ALL_TFRbs{4}.powspctrm(1,:,:,tbin1:tbin2),2)));    
        
        %absolute
        TFR_sessions(f,1,:,:) = squeeze(squeeze(mean(ALL_TFR{1}.powspctrm(1,:,:,tbin1:tbin2),2)));
        TFR_sessions(f,2,:,:) = squeeze(squeeze(mean(ALL_TFR{2}.powspctrm(1,:,:,tbin1:tbin2),2)));
        TFR_sessions(f,3,:,:) = squeeze(squeeze(mean(ALL_TFR{3}.powspctrm(1,:,:,tbin1:tbin2),2)));
        TFR_sessions(f,4,:,:) = squeeze(squeeze(mean(ALL_TFR{4}.powspctrm(1,:,:,tbin1:tbin2),2)));         

        %Pre-propofol
        TFR_prePropofol_bs_sessions(f,1,:,:) = squeeze(squeeze(mean(ALL_TFR{1}.powspctrm(1,:,:,tbin1PP:tbin2PP),2)));
        TFR_prePropofol_bs_sessions(f,2,:,:) = squeeze(squeeze(mean(ALL_TFR{2}.powspctrm(1,:,:,tbin1PP:tbin2PP),2)));
        TFR_prePropofol_bs_sessions(f,3,:,:) = squeeze(squeeze(mean(ALL_TFR{3}.powspctrm(1,:,:,tbin1PP:tbin2PP),2)));
        TFR_prePropofol_bs_sessions(f,4,:,:) = squeeze(squeeze(mean(ALL_TFR{4}.powspctrm(1,:,:,tbin1PP:tbin2PP),2)));        
        freqaxis = ALL_TFR{1}.freq;
    end
        
    drugstart_sessions(f) = time_axis_ROC(TFRbinstart);
    drugstart2_sessions(f) = time_axis_ROC(TFRbinstart2);
    drugstop_sessions(f) = time_axis_ROC(TFRbinstop);
    LOC_sessions(f) = time_axis_ROC(TFRbinLOC);
    ROC_sessions(f) = time_axis_ROC(TFRbinROC);
    
    
    
    clear TFR TFRbs ALL_TFR ALL_TFRbs sessionInfo
    
    f
    
end

% time_axis_ROC = time_axis_ROC(tbin1:tbin2);

time_axis_ROC = -30:1/6:9;

    xmarkers = 1:20:length(time_axis_ROC);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis_ROC(xmarkers(indx)));
        indx=indx+1;
    end
    
    
    %create frequency axis
    ymarkers = 1:5:length(freqaxis);
    ylabels = {};
    indx=1;
    for y = ymarkers
        ylabels{indx} = num2str(freqaxis(ymarkers(indx)), '%0.3g');
        indx=indx+1;
    end
    
%     sf1 = nearest(freqaxis, 0.1);
%     sf2 = nearest(freqaxis, 3);
%     
%     bs1 = nearest(time_axis_global,-60);
%     bs2 = nearest(time_axis_global,0);
%     
%     newTFRbs_sessions = squeeze((mean(TFRbs_sessions(:,:,sf1:sf2,:),3) - mean(mean(TFRbs_sessions(:,:,sf1:sf2,bs1:bs2),4),3)) ...
%           ./ mean(mean(TFRbs_sessions(:,:,sf1:sf2,bs1:bs2),4),3));
%       
% %       newTFRbs_sessions =  squeeze(mean(mean(TFRbs_sessions(:,:,sf1:sf2,bs1:bs2),4),3));
%     
%     figure; plot(time_axis_global, squeeze(nanmean(newTFRbs_sessions)))
%     xlim([-60 120])
    
    tmp = squeeze(nanmean(TFRbs_sessions(:,1,:,:)));
figure; imagesc(squeeze(nanmean(TFRbs_sessions(:,1,:,:))))

    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(time_axis_ROC, mean(drugstop_sessions)/60), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_ROC, 0), 'w', 'LineWidth', 2)
    caxis([-10 10])
    title('Manual dB')
%     xlim([nearest(time_axis_LOC, -60), nearest(time_axis_LOC, 120)])
    colorbar
    title('PFC')
    

figure; imagesc(squeeze(nanmean(TFRbs_sessions(:,2,:,:))))

    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(time_axis_ROC, mean(drugstop_sessions)/60), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_ROC, 0), 'w', 'LineWidth', 2)
    caxis([-10 10])
    title('Manual dB')
%     xlim([nearest(time_axis_ROC, -60), nearest(time_axis_ROC, 120)])
    colorbar
    title('8A')
    
        

figure; imagesc(squeeze(nanmean(TFRbs_sessions(:,3,:,:))))

    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(time_axis_ROC, mean(drugstop_sessions)/60), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_ROC, 0), 'w', 'LineWidth', 2)
    caxis([-10 10])
    title('Manual dB')
%     xlim([nearest(time_axis_ROC, -60), nearest(time_axis_ROC, 120)])
    colorbar
    title('PPC')
    
figure; imagesc(squeeze(nanmean(TFRbs_sessions(:,4,:,:))))

    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(time_axis_ROC, mean(drugstop_sessions)/60), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_ROC, 0), 'w', 'LineWidth', 2)
    caxis([-10 10])
    title('Manual dB')
%     xlim([nearest(time_axis_ROC, -60), nearest(time_axis_ROC, 120)])
    colorbar
    title('STG')
    
    
    
    return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Statistics
%%
% Channel level randomization

arealist = {'PFC', '8A', 'PPC', 'STG'};    

numfreq=length(freqaxis);
time_axis_LOC = -20:1/6:-15;

numtime = length(time_axis_ROC);
bs1 = 1
bs2 = nearest(time_axis_LOC,-15);
first_level_crit = 0.01;

% All_sigclust = zeros( numfreq, numtime,5);
% All_empirical_effect = zeros(  numfreq, numtime,5);
nrand = 1000;
rand_clust_size = zeros(4,nrand);
for a =1:4
    
    numsess = size(TFRbs_sessions,1);
    All_empirical_effect(:,:,a) = squeeze(mean(TFRbs_sessions(a,:,:,:)));
    %perform randomization
    for r = 1:nrand
        TFR_rand = squeeze(cat(4, TFR_prePropofol_bs_sessions, TFRbs_sessions)); %num sess x numfreq x numtime x numarea
        numtime = size(TFR_rand,4);
        num_bs = length(bs1:bs2);
        num_post = length(bs2+1:numtime);

        for aa = 1:numsess
            %randomize frequency and pre vs. post laser power
            baseline_bins = repmat(1:num_bs, [1 10]);
            baseline_bins = baseline_bins(randperm(length(baseline_bins)));
            baseline_bins = baseline_bins(1:num_post);
            drug_bins = bs2+1:numtime;
            allbins = [baseline_bins drug_bins];
            allbins = allbins(randperm(length(allbins)));
            rand_bins = allbins(1:numtime);
            
            TFR_rand(aa,a,:,:) = TFR_rand(aa,a,:,rand_bins);  
        end

        first_level_stat = zeros(numfreq,numtime);
        for tbin = 1:numtime
            [h p] = ttest2(squeeze(mean(TFR_rand(:,a,:,bs1:bs2),4)), squeeze(TFR_rand(:,a,:,tbin)));
            first_level_stat(:,tbin) = p;
%             tbin
        end
        
        random_effect =  squeeze(first_level_stat(:,:)) < first_level_crit;
        random_mat_pass = zeros(numfreq, numtime);
        random_mat_pass(find(random_effect==1)) = 1;
        
        [L,num] = spm_bwlabel(random_mat_pass,6);
        Lpass = zeros(size(L));
        clustsize = zeros(num,1);

        for clus = 1:num
            clustsize(clus) = length(find(L==clus));
        end
        
        rand_clust_size(a,r) = max(clustsize);
        
        
        r
    end
    
end

cd \\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis
save TFR_ROC_locked_stats rand_clust_size TFR_prePropofol_bs_sessions TFRbs_sessions time_axis_ROC freqaxis first_level_crit nrand drugstop_sessions

%%
cd \\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis
load TFR_ROC_locked_stats
numtime = length(time_axis_ROC);
arealist = {'PFC', '8A', 'PPC', 'STG'};    

numfreq=length(freqaxis);
time_axis_LOC = -20:1/6:-15;
bs1 = 1
bs2 = nearest(time_axis_LOC,-15);

corrected_alphaval = 0.01;    
for a=1:4
    
    sorted_clusters = sort(rand_clust_size(a,:), 'descend');
    alphabin = floor(corrected_alphaval.*nrand);
    threshold = sorted_clusters(alphabin);

    first_level_stat = zeros(numfreq,numtime);
    %get empirical results
    numtime = size(TFRbs_sessions,4);
    for tbin = 1:numtime
        [h p] = ttest2(squeeze(mean(TFR_prePropofol_bs_sessions(:,a,:,bs1:bs2),4)), squeeze(TFRbs_sessions(:,a,:,tbin)));
        first_level_stat(:,tbin) = p;
        %             tbin
    end
    
    empirical_effect =  squeeze(first_level_stat(:,:)) < first_level_crit;
    empirical_mat_pass = zeros(numfreq, numtime);
    empirical_mat_pass(find(empirical_effect==1)) = 1;
    
    [L,num] = spm_bwlabel(empirical_mat_pass,6);
    Lpass = 0.5.*ones(size(L));
    clustsize = zeros(num,1);
    sigclussize = [];
    for clus = 1:num
        clustsize(clus) = length(find(L==clus));
        if clustsize(clus) > threshold
            Lpass(L==clus) = 1;
            %                 hold on; plot([clustsize(clus) clustsize(clus)], [0 1000], 'r--')
        end
    end
    All_sigclust(:,:,a) = Lpass;
    
    %create time axis
    xmarkers = 1:20:length(time_axis_ROC);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis_ROC(xmarkers(indx)));
        indx=indx+1;
    end
    
    %create frequency axis
    ymarkers = 1:5:length(freqaxis);
    ylabels = {};
    indx=1;
    for y = ymarkers
        ylabels{indx} = num2str(freqaxis(ymarkers(indx)), '%0.3g');
        indx=indx+1;
    end    
    
    
    figure; imagesc(squeeze(mean(TFRbs_sessions(:,a,:,:,:))));
    set(gca, 'ydir', 'normal')
    caxis([-10 10])
    colorbar;
    
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since ROC (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')    
    alpha(Lpass)
    hold on;
    contourLines = zeros(size(Lpass));
    contourLines(find(Lpass==1)) = 1;
    contour(contourLines,1, 'k')
    xline(nearest(time_axis_ROC, 0), 'k--')
    title(arealist{a})
    xlim([nearest(time_axis_ROC, -20) nearest(time_axis_ROC, 10)])
%     xline(nearest(time_axis_ROC, mean(drugstart_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_ROC, mean(drugstop_sessions)/60), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_ROC, mean(drugstop_sessions)/60)+(std(drugstop_sessions)/60), 'k', 'LineWidth', 1)
    xline(nearest(time_axis_ROC, mean(drugstop_sessions)/60)-(std(drugstop_sessions)/60), 'k', 'LineWidth', 1)
%     xline(nearest(time_axis_ROC, mean(LOC_sessions)), 'w', 'LineWidth', 2)
%     xline(nearest(time_axis_ROC, mean(ROC_sessions)), 'w', 'LineWidth', 2)    
    
end
    
    
    
    