

cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')


filestoload = ls('*_TFR_spikes*')

spike_rates_sessions = {};
time_axis_global = -80:10/60:180;
TFRbs_sessions = NaN(size(filestoload,1), 4,50,length(-20*60:10:30*60));
TFR_sessions = NaN(size(filestoload,1), 4,50,length(-20*60:10:30*60));
nunits_global = ones(4,1);
drugstart_sessions = zeros(size(filestoload,1),1);
drugstart2_sessions = zeros(size(filestoload,1),1);
drugstop_sessions = zeros(size(filestoload,1),1);
LOC_sessions = zeros(size(filestoload,1),1);
ROC_sessions = zeros(size(filestoload,1),1);
for f = 1:size(filestoload,1)

    load(filestoload(f,:)); %, 'sessionInfo', 'spike_rates', 'time_axis_baseline')
    
    if isempty(sessionInfo.eyesClose)
        error('deal with no eye info')
    else
        time_axis_LOC = time_axis_baseline-((sessionInfo.eyesClose(1)-sessionInfo.drugStart(1))/60);
    end
   
    
    tbin1 = nearest(time_axis_LOC, -20);
    tbin2 = nearest(time_axis_LOC, 30);
    
    for a = 1:4
        nunits = size(spike_rates{a},1);
%         spike_rates_sessions{a}(nunits_global(a):nunits_global(a)+nunits-1, tbin1:tbin2) = NaN;
        spike_rates_sessions{a}(nunits_global(a):nunits_global(a)+nunits-1, :) = spike_rates{a}(:,tbin1:tbin2);
    
        nunits_global(a) = nunits_global(a) + nunits;
    end
    
    
    if f<=11 
        
    TFRbinstart = nearest(TFR.time, sessionInfo.drugStart(1));
    TFRbinstart2 = nearest(TFR.time, sessionInfo.drugStart(2));
    TFRbinstop = nearest(TFR.time, sessionInfo.drugStart(end)+(sessionInfo.drugDuration(end)*60));
    TFRbinLOC = nearest(TFR.time, sessionInfo.eyesClose(1));
    TFRbinROC = nearest(TFR.time, sessionInfo.eyesOpen(1));        
        
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
    else        

    TFRbinstart = nearest(ALL_TFR{1}.time, sessionInfo.drugStart(1));
    TFRbinstart2 = nearest(ALL_TFR{1}.time, sessionInfo.drugStart(2));
    TFRbinstop = nearest(ALL_TFR{1}.time, sessionInfo.drugStart(end)+(sessionInfo.drugDuration(end)*60));
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
        
        freqaxis = ALL_TFR{1}.freq;
    end
        
    drugstart_sessions(f) = time_axis_LOC(TFRbinstart);
    drugstart2_sessions(f) = time_axis_LOC(TFRbinstart2);
    drugstop_sessions(f) = time_axis_LOC(TFRbinstop);
    LOC_sessions(f) = time_axis_LOC(TFRbinLOC);
    ROC_sessions(f) = time_axis_LOC(TFRbinROC);
    
    
    
    clear TFR TFRbs ALL_TFR ALL_TFRbs sessionInfo
    
    f
    
end

time_axis_LOC = time_axis_LOC(tbin1:tbin2);

    xmarkers = 18:25:length(time_axis_LOC);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis_LOC(xmarkers(indx)));
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
    xline(nearest(time_axis_LOC, mean(drugstart_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(drugstop_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(LOC_sessions)), 'w', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(ROC_sessions)), 'w', 'LineWidth', 2)
%     caxis([-10 10])
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
    xline(nearest(time_axis_LOC, mean(drugstart_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(drugstop_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(LOC_sessions)), 'w', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(ROC_sessions)), 'w', 'LineWidth', 2)
    caxis([-10 10])
    title('Manual dB')
    xlim([nearest(time_axis_LOC, -60), nearest(time_axis_LOC, 120)])
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
    xline(nearest(time_axis_LOC, mean(drugstart_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(drugstop_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(LOC_sessions)), 'w', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(ROC_sessions)), 'w', 'LineWidth', 2)
    caxis([-10 10])
    title('Manual dB')
    xlim([nearest(time_axis_LOC, -60), nearest(time_axis_LOC, 120)])
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
    xline(nearest(time_axis_LOC, mean(drugstart_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(drugstop_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(LOC_sessions)), 'w', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(ROC_sessions)), 'w', 'LineWidth', 2)
    caxis([-10 10])
    title('Manual dB')
    xlim([nearest(time_axis_LOC, -60), nearest(time_axis_LOC, 120)])
    colorbar
    title('STG')
    
    return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Statistics
%%
% Channel level randomization

arealist = {'PFC', '8A', 'PPC', 'STG'};    

numfreq=length(freqaxis);

numtime = length(time_axis_LOC);
bs1 = 1;
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
        TFR_rand = squeeze(TFRbs_sessions(:,:,:,:)); %num sess x numfreq x numtime x numarea
        
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
save TFR_LOC_locked_stats rand_clust_size TFRbs_sessions time_axis_LOC freqaxis first_level_crit nrand drugstart_sessions


%%
cd \\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis
load TFR_LOC_locked_stats
numtime = length(time_axis_LOC);
arealist = {'PFC', '8A', 'PPC', 'STG'};    
bs1 = 1;
bs2 = nearest(time_axis_LOC,-15);
numfreq=length(freqaxis);
corrected_alphaval = 0.01;    
for a=1:4
    
    sorted_clusters = sort(rand_clust_size(a,:), 'descend');
    alphabin = floor(corrected_alphaval.*nrand);
    threshold = sorted_clusters(alphabin);

    first_level_stat = zeros(numfreq,numtime);
    %get empirical results
    for tbin = 1:numtime
        [h p] = ttest2(squeeze(mean(TFRbs_sessions(:,a,:,bs1:bs2),4)), squeeze(TFRbs_sessions(:,a,:,tbin)));
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
    xmarkers = 1:20:length(time_axis_LOC);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis_LOC(xmarkers(indx)));
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
    xlabel('Time since LOC (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')    
    
    alpha(Lpass)
    hold on;
    contourLines = zeros(size(Lpass));
    contourLines(find(Lpass==1)) = 1;
    contour(contourLines,1, 'k')
    xline(nearest(time_axis_LOC, 0), 'k--')
    title(arealist{a})
    xlim([nearest(time_axis_LOC, -15) nearest(time_axis_LOC, 20)])
    xline(nearest(time_axis_LOC, mean(drugstart_sessions)), 'k', 'LineWidth', 2)
    xline(nearest(time_axis_LOC, mean(drugstart_sessions))+(std(drugstart_sessions)), 'k', 'LineWidth', 1)
    xline(nearest(time_axis_LOC, mean(drugstart_sessions))-(std(drugstart_sessions)), 'k', 'LineWidth', 1)
%     xline(nearest(time_axis_LOC, mean(drugstop_sessions)), 'k', 'LineWidth', 2)
%     xline(nearest(time_axis_LOC, mean(LOC_sessions)), 'w', 'LineWidth', 2)
%     xline(nearest(time_axis_LOC, mean(ROC_sessions)), 'w', 'LineWidth', 2)    
    
end
    
%%%%%%%%%%%%%%%

%%

%Compare spetra between areas

numtime = length(time_axis_LOC);
arealist = {'PFC', '8A', 'PPC', 'STG'};    
bs1 = 1;
bs2 = nearest(time_axis_LOC,-15);

figure; 
loglog(freqaxis, squeeze(mean(mean(TFR_sessions(:,1,:,bs1:bs2),4),1)), 'b'); hold on
loglog(freqaxis, squeeze(mean(mean(TFR_sessions(:,2,:,bs1:bs2),4),1)), 'r')
loglog(freqaxis, squeeze(mean(mean(TFR_sessions(:,3,:,bs1:bs2),4),1)), 'g'); hold on
loglog(freqaxis, squeeze(mean(mean(TFR_sessions(:,4,:,bs1:bs2),4),1)), 'm')
bs1 = nearest(time_axis_LOC,10);
bs2 = nearest(time_axis_LOC,30);

loglog(freqaxis, squeeze(mean(mean(TFR_sessions(:,1,:,bs1:bs2),4),1)), 'b--')
loglog(freqaxis, squeeze(mean(mean(TFR_sessions(:,2,:,bs1:bs2),4),1)), 'r--')
loglog(freqaxis, squeeze(mean(mean(TFR_sessions(:,3,:,bs1:bs2),4),1)), 'g--'); hold on
loglog(freqaxis, squeeze(mean(mean(TFR_sessions(:,4,:,bs1:bs2),4),1)), 'm--')


nsess = size(TFRbs_sessions,1);
sem1 = squeeze(std(mean(TFRbs_sessions(:,1,:,bs1:bs2),4),1))./sqrt(nsess);
sem2 = squeeze(std(mean(TFRbs_sessions(:,2,:,bs1:bs2),4),1))./sqrt(nsess);
sem3 = squeeze(std(mean(TFRbs_sessions(:,3,:,bs1:bs2),4),1))./sqrt(nsess);
sem4 = squeeze(std(mean(TFRbs_sessions(:,4,:,bs1:bs2),4),1))./sqrt(nsess);
figure
shadedErrorBar(freqaxis, squeeze(mean(mean(TFRbs_sessions(:,1,:,bs1:bs2),4),1)), sem1, 'b'); hold on
[PKS,LOCS,W] = findpeaks(squeeze(mean(mean(TFRbs_sessions(:,1,:,bs1:bs2),4),1)), freqaxis)
xline(LOCS(1), 'b--')
shadedErrorBar(freqaxis, squeeze(mean(mean(TFRbs_sessions(:,2,:,bs1:bs2),4),1)), sem2, 'r');
[PKS,LOCS,W] = findpeaks(squeeze(mean(mean(TFRbs_sessions(:,2,:,bs1:bs2),4),1)), freqaxis)
xline(LOCS(1), 'r--')
xline(LOCS(2), 'r--')
shadedErrorBar(freqaxis, squeeze(mean(mean(TFRbs_sessions(:,3,:,bs1:bs2),4),1)), sem3, 'g'); hold on
[PKS,LOCS,W] = findpeaks(squeeze(mean(mean(TFRbs_sessions(:,3,:,bs1:bs2),4),1)), freqaxis)
xline(LOCS(1), 'g--')
shadedErrorBar(freqaxis, squeeze(mean(mean(TFRbs_sessions(:,4,:,bs1:bs2),4),1)), sem4, 'm');
[PKS,LOCS,W] = findpeaks(squeeze(mean(mean(TFRbs_sessions(:,4,:,bs1:bs2),4),1)), freqaxis)
xline(LOCS(1), 'm--')
yline(0, 'k--')
set(gca, 'xscale', 'log')
xlabel('Frequency (Hz)')
ylabel('dB Change in power')
xlim([freqaxis(1) 200])


%%
%test for significance of peak differences


bs1 = nearest(time_axis_LOC,10);
bs2 = nearest(time_axis_LOC,30);

nsess = size(TFRbs_sessions,1);
peaks = zeros(nsess, 4);
for f = 1:nsess
    %PFC
    spectrum1 = squeeze(mean(TFRbs_sessions(f,1,:,bs1:bs2),4));
    spectrum1 = smooth(spectrum1, 5);
    [PKS,LOCS1,W] = findpeaks(spectrum1, freqaxis, 'NPeaks', 1)
    %8A
    spectrum2 = squeeze(mean(TFRbs_sessions(f,2,:,bs1:bs2),4));
    spectrum2 = smooth(spectrum2, 5);
    [PKS,LOCS2,W] = findpeaks(spectrum2, freqaxis, 'NPeaks', 1)
    %PPC
    spectrum3 = squeeze(mean(TFRbs_sessions(f,3,:,bs1:bs2),4));
    spectrum3 = smooth(spectrum3, 5);
    [PKS,LOCS3,W] = findpeaks(spectrum3, freqaxis, 'NPeaks', 1)
    %STG
    spectrum4 = squeeze(mean(TFRbs_sessions(f,4,:,bs1:bs2),4));
    spectrum4 = smooth(spectrum4, 5);
    [PKS,LOCS4,W] = findpeaks(spectrum4, freqaxis, 'NPeaks', 1)    
%     rmi = find(PKS<2)
%     PKS(rmi) = [];
%     LOCS(rmi) = [];
%     W(rmi) = [];
    peaks(f,1) = LOCS1(1);
    peaks(f,2) = LOCS2(1);
    peaks(f,3) = LOCS3(1);
    peaks(f,4) = LOCS4(1);
    figure; 
    plot(freqaxis, spectrum1, 'b'); hold on;
    plot(freqaxis, spectrum2, 'r'); hold on;
    plot(freqaxis, spectrum3, 'g'); hold on;
    plot(freqaxis, spectrum4, 'm'); hold on;
    set(gca, 'xscale', 'log')
%     for p = 1:length(PKS)
    xline(LOCS1(1), 'b--')
    xline(LOCS2(1), 'r--')
    xline(LOCS3(1), 'g--')
    xline(LOCS4(1), 'm--')
%     end
end

[p h] = signtest(peaks(:,1)-peaks(:,3))
[p h] = signtest(peaks(:,1)-peaks(:,4))
[p h] = signtest(peaks(:,2)-peaks(:,3))
[p h] = signtest(peaks(:,2)-peaks(:,4))

%all sig at p = 1.9073e-06

%%
%Anteriorization
bs1 = 1;
bs2 = nearest(time_axis_LOC,-15);
frontal_vs_posterior_pow_awake = squeeze(mean(mean(TFR_sessions(:,1:2,:,bs1:bs2),4),2)) ./ ...
    squeeze(mean(mean(TFR_sessions(:,3:4,:,bs1:bs2),4),2));


bs1 = nearest(time_axis_LOC,10);
bs2 = nearest(time_axis_LOC,30);
frontal_vs_posterior_pow_unconsc = squeeze(mean(mean(TFR_sessions(:,1:2,:,bs1:bs2),4),2)) ./ ...
    squeeze(mean(mean(TFR_sessions(:,3:4,:,bs1:bs2),4),2));


sem1 = std(frontal_vs_posterior_pow_awake)./sqrt(nsess);
sem2 = std(frontal_vs_posterior_pow_unconsc)./sqrt(nsess);

figure
shadedErrorBar(freqaxis, mean(frontal_vs_posterior_pow_awake), sem1, 'b'); hold on
shadedErrorBar(freqaxis, mean(frontal_vs_posterior_pow_unconsc), sem2, 'k');

yline(1, 'k--')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlim([0.178 200])

%%

b1 = nearest(freqaxis, 12);
b2 = nearest(freqaxis, 25);

bs1 = nearest(time_axis_LOC,10);
bs2 = nearest(time_axis_LOC,30);

pow_post = squeeze(mean(mean(TFR_sessions(:,1:4,b1:b2,bs1:bs2),4),3));
bs1 = nearest(time_axis_LOC,-20);
bs2 = nearest(time_axis_LOC,-15);

pow_pre = squeeze(mean(mean(TFR_sessions(:,1:4,b1:b2,bs1:bs2),4),3));
figure; 
subplot(2,1,1);
bar(mean(pow_pre));
hold on;
errorbar(1:4, mean(pow_pre), std(pow_pre)./nsess, 'linestyle', 'none')
subplot(2,1,2);
bar(mean(pow_post))
hold on;
errorbar(1:4, mean(pow_post), std(pow_post)./nsess, 'linestyle', 'none')

[p h] = signtest(pow_pre(:,1)-pow_pre(:,3))
[p h] = signtest(pow_pre(:,1)-pow_pre(:,4))
[p h] = signtest(pow_pre(:,2)-pow_pre(:,3))
[p h] = signtest(pow_pre(:,2)-pow_pre(:,4))

[p h] = signtest(pow_post(:,1)-pow_post(:,3))
[p h] = signtest(pow_post(:,1)-pow_post(:,4))
[p h] = signtest(pow_post(:,2)-pow_post(:,3))
[p h] = signtest(pow_post(:,2)-pow_post(:,4))

