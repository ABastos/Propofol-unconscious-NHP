
cd \\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolWakeUp

% filestoload = ls('\\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolWakeUp\*.mat');
filestoload = ls('\\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolWakeUp_sorted\*.mat');

% filestoload(end,:) = [];

for s = 17:20 % % 1:size(filestoload,1)
% for s = 11:size(filestoload,1)
    cd('\\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolWakeUp_sorted\')
    load(deblank(filestoload(s,:)), 'spikeWaves', 'trialInfo', 'sessionInfo', 'electrodeInfo', 'spikeChnlInfo', 'spikeTimes', 'spikeTimesSchema', 'unitInfo')
        
    
    isStandardAllAmps = trialInfo.dbs_isWakeUpTest & ...      % Is it a wake-up test?
        ismember(trialInfo.dbs_testType,{'WUT7','WUT7ab','WUT7abU'}) & ... % WUT7/ab family of tests
        ... (trialInfo.dbs_amplitude > 0.5) & ... % Amplitude > 0.5 mAmp
        (trialInfo.dbs_frequency == 180) & ...% Frequency = 180 hz
        (trialInfo.dbs_duration > 25) & ...   % 20 s < Duration < 35 s
        (trialInfo.dbs_duration < 35) & ...
        (trialInfo.dbs_previousITI > 120) ; %previous stim was greater than 2 minutes ago
    
    standardTrialIndx = find(isStandardAllAmps);
    
    numtr = length(standardTrialIndx);
    
    
    %Get spike rate modulation
    
    %find units with at least 0.1 Hz firing rate over session
%     units = length(find(unitInfo.numSpikes > 1000));
    
    %find un-sorted units in each array
    PFCindx2 = find(strcmp(unitInfo.array, 'PFC') & unitInfo.unit == 0);
    FEFindx2 = find(strcmp(unitInfo.array, 'FEF') & unitInfo.unit == 0);
    PPCindx2 = find(strcmp(unitInfo.array, 'PPC') & unitInfo.unit == 0);
    STGindx2 = find(strcmp(unitInfo.array, 'STG') & unitInfo.unit == 0);
    allunits2 = [PFCindx2; FEFindx2; PPCindx2; STGindx2];
    
    %only sorted spikes
    PFCindx = find(strcmp(unitInfo.array, 'PFC') & unitInfo.unit > 0);
    FEFindx = find(strcmp(unitInfo.array, 'FEF') & unitInfo.unit > 0);
    PPCindx = find(strcmp(unitInfo.array, 'PPC') & unitInfo.unit > 0);
    STGindx = find(strcmp(unitInfo.array, 'STG') & unitInfo.unit > 0);    
    
    allunits = [PFCindx; FEFindx; PPCindx; STGindx];
    
    %find the max time bin in the session:
    maxtime = zeros(length(spikeTimes),1);
    for c = 1:length(spikeTimes)
       maxtime(c) = floor(max(spikeTimes{c})*1000);
    end
    ntime = max(maxtime);
    
    spike_waveforms_pre = zeros(length(allunits), 48);
    spike_waveforms_post = zeros(length(allunits), 48);
    spike_waveforms_stim = zeros(length(allunits), 48);
    spike_waveforms_NOstim = zeros(length(allunits), 48);
    for c = 1:length(allunits)
        %get all spikes that occurred prior to propofol administration
        indx_pre = find(spikeTimes{allunits(c)} < sessionInfo.drugStart(1));
        indx_post = find(spikeTimes{allunits(c)} > sessionInfo.drugEnd(end));
        indx_stim = [];
        for tr = 1:numtr
            indx_stim_tmp = find(spikeTimes{allunits(c)} >= trialInfo.dbs_stimOn(standardTrialIndx(tr)) & ...
                             spikeTimes{allunits(c)} <= trialInfo.dbs_stimOff(standardTrialIndx(tr)));
            indx_stim = [indx_stim indx_stim_tmp];
        end
        indx_NOstim = setdiff(1:length(spikeTimes{allunits(c)}), [indx_pre'; indx_post'; indx_stim']);
        
        spike_waveforms_pre(c,:) = mean(spikeWaves{allunits(c)}(:,indx_pre),2);
        spike_waveforms_stim(c,:) = mean(spikeWaves{allunits(c)}(:,indx_stim),2);
        spike_waveforms_NOstim(c,:) = mean(spikeWaves{allunits(c)}(:,indx_NOstim),2);
        spike_waveforms_post(c,:) = mean(spikeWaves{allunits(c)}(:,indx_post),2);
    end    
    
    corr_values = zeros(length(allunits),3);
%     figure;
    for c = 1:length(allunits)
        [r3 p] = corr(spike_waveforms_pre(c,:)', spike_waveforms_post(c,:)', 'type', 'Pearson');
        [r2 p] = corr(spike_waveforms_pre(c,:)', spike_waveforms_stim(c,:)', 'type', 'Pearson');
        [r1 p] = corr(spike_waveforms_pre(c,:)', spike_waveforms_NOstim(c,:)', 'type', 'Pearson');
        corr_values(c,1) = r1;
        corr_values(c,2) = r2;
        corr_values(c,3) = r3;
%         plot(1:48, spike_waveforms_pre(c,:), 'b'); hold on;
%         plot(1:48, spike_waveforms_post(c,:), 'k');
%         plot(1:48, spike_waveforms_stim(c,:), 'r');
%         plot(1:48, spike_waveforms_NOstim(c,:), 'r--', 'lineWidth', 2);
%         title(sprintf('Unit: %i, pre vs. post = %0.5g, pre vs. stim = %0.5g, pre vs. NOstim = %0.5g', c, r3,r2,r1))
%         pause
%         clf
    end
    
    figure;
    subplot(3,1,1);
    hist(corr_values(:,1),40); hold on; xlim([0.7 1]); xline(0.99, 'r'); title('Pre vs Post'); 
    subplot(3,1,2);
    hist(corr_values(:,2),40); hold on; xlim([0.7 1]); xline(0.99, 'r'); title('Pre vs Stim'); 
    subplot(3,1,3);
    hist(corr_values(:,3),40); hold on; xlim([0.7 1]); xline(0.99, 'r'); title('Pre vs NoStim');
    
    keep_units = find(corr_values(:,1) > 0.99 & corr_values(:,2) > 0.99 & corr_values(:,3) > 0.99);
    
    disp(['keeping ' int2str(length(keep_units)) ' out of ' int2str(length(allunits)) ' units'])
    kept_vs_all = [length(keep_units) length(allunits)];
    
    allunits = allunits(keep_units);

    FT_data = [];
    FT_data.trial{1} = zeros(length(allunits), ntime);
    FT_data.label = {};
    for c = 1:length(allunits)
        FT_data.label{c} = [unitInfo.array{allunits(c)} '-' int2str(c)];
        timeIndx = floor(spikeTimes{allunits(c)}.*1000);
        rmindx = find(timeIndx>ntime);
        timeIndx(rmindx) = [];
        spikeTrain = timeIndx;
        FT_data.trial{1}(c,spikeTrain) = 1;
    end    

    %get data into DBS trial format
    data =[];
    data.trialinfo = zeros(numtr,4);
    for tr = 1:numtr
        data.time{tr} = -30:0.001:120;
        data.trial{tr} = zeros(length(allunits), length(data.time{1})); %num channel x num time
        tr_start = trialInfo.dbs_stimOn(standardTrialIndx(tr))*1000;
        data.trial{tr}(:,:) = squeeze(FT_data.trial{1}(:,tr_start-(30*1000):tr_start+(120*1000)));
        data.trialinfo(tr,1) = standardTrialIndx(tr);
        data.trialinfo(tr,2) = trialInfo.dbs_amplitude(standardTrialIndx(tr));
        data.trialinfo(tr,3) = trialInfo.dbs_wakeScore(standardTrialIndx(tr));
        data.trialinfo(tr,4) = trialInfo.dbs_duration(standardTrialIndx(tr));
    end
    data.label = FT_data.label;
    data.fsample = 1000;
    if strcmp(sessionInfo.subject, 'Mary')
        data.trialinfo(:,5) = 2.*ones(numtr,1);
    else
        data.trialinfo(:,5) = ones(numtr,1);
    end
    
%     %now for the unsorted spikes
%     FT_data2 = [];
%     FT_data2.trial{1} = zeros(length(allunits2), ntime);
%     FT_data2.label = {};
%     for c = 1:length(allunits2)
%         FT_data2.label{c} = [unitInfo.array{allunits2(c)} '-' int2str(c)];
%         timeIndx = floor(spikeTimes{allunits2(c)}.*1000);
%         rmindx = find(timeIndx>ntime);
%         timeIndx(rmindx) = [];
%         spikeTrain = timeIndx;
%         FT_data2.trial{1}(c,spikeTrain) = 1;
%     end    
% 
%     %get data into DBS trial format
%     data2 =[];
%     data2.trialinfo = zeros(numtr,4);
%     for tr = 1:numtr
%         data2.time{tr} = -30:0.001:120;
%         data2.trial{tr} = zeros(length(allunits2), length(data.time{1})); %num channel x num time
%         tr_start = trialInfo.dbs_stimOn(standardTrialIndx(tr))*1000;
%         data2.trial{tr}(:,:) = squeeze(FT_data2.trial{1}(:,tr_start-(30*1000):tr_start+(120*1000)));
%         data2.trialinfo(tr,1) = standardTrialIndx(tr);
%         data2.trialinfo(tr,2) = trialInfo.dbs_amplitude(standardTrialIndx(tr));
%         data2.trialinfo(tr,3) = trialInfo.dbs_wakeScore(standardTrialIndx(tr));
%         data2.trialinfo(tr,4) = trialInfo.dbs_duration(standardTrialIndx(tr));
%     end
%     data2.label = FT_data.label;
%     data2.fsample = 1000;
%     if strcmp(sessionInfo.subject, 'Mary')
%         data2.trialinfo(:,5) = 2.*ones(numtr,1);
%     else
%         data2.trialinfo(:,5) = ones(numtr,1);
%     end
    
    
    
    newTrialInfo = data.trialinfo;

%     cfg = [];
%     cfg.viewmode = 'vertical';
% %     cfg.channel = [1 size(electrodeInfo.electrode,1)+1:length(data.label)];
% %     cfg.channel = [1 65 126 190];
%     % cfg.chanscale = stdchan([1 size(electrodeInfo.electrode,1)+1:length(data.label)]);
%     
%     ft_databrowser(cfg, data)    
    
    
    PFCindx = [];
    FEFindx = [];
    PPCindx = [];
    STGindx = [];
    for c = 1:length(data.label)
        if strcmp(data.label{c}(1:3), 'PPC')
            PPCindx = [PPCindx; c];
        elseif strcmp(data.label{c}(1:3), 'FEF')
            FEFindx = [FEFindx; c];
        elseif strcmp(data.label{c}(1:3), 'STG')
            STGindx = [STGindx; c];
        elseif strcmp(data.label{c}(1:3), 'PFC')
            PFCindx = [PFCindx; c];
        end
    end
    
    %Get psth
    time_axis = -30:1:120;
    winsize = 2;
    %trial by unit by time
    
    area_labels = {'PFC', 'FEF', 'PPC', 'STG'};
    
    spike_rates = cell(4,1);
    spike_rates{1} = zeros(length(data.trial), length(PFCindx),length(time_axis));
    spike_rates{2} = zeros(length(data.trial), length(FEFindx),length(time_axis));
    spike_rates{3} = zeros(length(data.trial), length(PPCindx),length(time_axis));
    spike_rates{4} = zeros(length(data.trial), length(STGindx),length(time_axis));
    
    spike_rates_bs = cell(4,1);
    spike_rates_bs{1} = zeros(length(data.trial), length(PFCindx),length(time_axis));
    spike_rates_bs{2} = zeros(length(data.trial), length(FEFindx),length(time_axis));
    spike_rates_bs{3} = zeros(length(data.trial), length(PPCindx),length(time_axis));
    spike_rates_bs{4} = zeros(length(data.trial), length(STGindx),length(time_axis));
    
    for tr = 1:length(data.trial)
        for c = 1:length(PFCindx)
            for tbin = 1:length(time_axis)
                
                tbin1 = nearest(data.time{1}, time_axis(tbin)-winsize/2);
                tbin2 = nearest(data.time{1}, time_axis(tbin)+winsize/2);
                time1 = data.time{1}(tbin1);
                time2 = data.time{1}(tbin2);
                spike_rates{1}(tr,c,tbin) = sum(data.trial{tr}(PFCindx(c),tbin1:tbin2))/(time2-time1);
            end
        end
        bs1 = nearest(time_axis, -30);
        bs2 = nearest(time_axis, -2);
        
        spike_rates_bs{1}(tr,:,:) = spike_rates{1}(tr,:,:) - mean(spike_rates{1}(tr,:,bs1:bs2),3);
        tr
    end
    for tr = 1:length(data.trial)
        for c = 1:length(FEFindx)
            for tbin = 1:length(time_axis)
                
                tbin1 = nearest(data.time{1}, time_axis(tbin)-winsize/2);
                tbin2 = nearest(data.time{1}, time_axis(tbin)+winsize/2);
                time1 = data.time{1}(tbin1);
                time2 = data.time{1}(tbin2);
                spike_rates{2}(tr,c,tbin) = sum(data.trial{tr}(FEFindx(c),tbin1:tbin2))/(time2-time1);
            end
        end
        bs1 = nearest(time_axis, -30);
        bs2 = nearest(time_axis, -2);
        
        spike_rates_bs{2}(tr,:,:) = spike_rates{2}(tr,:,:) - mean(spike_rates{2}(tr,:,bs1:bs2),3);
        tr
    end
    for tr = 1:length(data.trial)
        for c = 1:length(PPCindx)
            for tbin = 1:length(time_axis)
                
                tbin1 = nearest(data.time{1}, time_axis(tbin)-winsize/2);
                tbin2 = nearest(data.time{1}, time_axis(tbin)+winsize/2);
                time1 = data.time{1}(tbin1);
                time2 = data.time{1}(tbin2);
                spike_rates{3}(tr,c,tbin) = sum(data.trial{tr}(PPCindx(c),tbin1:tbin2))/(time2-time1);
            end
        end
        bs1 = nearest(time_axis, -30);
        bs2 = nearest(time_axis, -2);
        
        spike_rates_bs{3}(tr,:,:) = spike_rates{3}(tr,:,:) - mean(spike_rates{3}(tr,:,bs1:bs2),3);
        tr
    end
    for tr = 1:length(data.trial)
        for c = 1:length(STGindx)
            for tbin = 1:length(time_axis)
                
                tbin1 = nearest(data.time{1}, time_axis(tbin)-winsize/2);
                tbin2 = nearest(data.time{1}, time_axis(tbin)+winsize/2);
                time1 = data.time{1}(tbin1);
                time2 = data.time{1}(tbin2);
                spike_rates{4}(tr,c,tbin) = sum(data.trial{tr}(STGindx(c),tbin1:tbin2))/(time2-time1);
                
            end
        end
        bs1 = nearest(time_axis, -30);
        bs2 = nearest(time_axis, -2);
        
        spike_rates_bs{4}(tr,:,:) = spike_rates{4}(tr,:,:) - mean(spike_rates{4}(tr,:,bs1:bs2),3);
        tr
    end
    
%     tbin1 = nearest(time_axis, -1);
%     tbin2 = nearest(time_axis, 30);
%     
%     spike_rates{1}(:,:,tbin1:tbin2) = NaN;
%     spike_rates{2}(:,:,tbin1:tbin2) = NaN;
%     spike_rates{3}(:,:,tbin1:tbin2) = NaN;
%     spike_rates{4}(:,:,tbin1:tbin2) = NaN;
%     
%     spike_rates_bs{1}(:,:,tbin1:tbin2) = NaN;
%     spike_rates_bs{2}(:,:,tbin1:tbin2) = NaN;
%     spike_rates_bs{3}(:,:,tbin1:tbin2) = NaN;
%     spike_rates_bs{4}(:,:,tbin1:tbin2) = NaN;
    

    if strcmp(sessionInfo.subject, 'Mary')
        %mary
        low_current = [0.5 1.0];
        high_current = [1.5 2.5];
%         tmp_trialinfo(:,5) = 2.*ones(size(ALL_TFRbs{1}.trialinfo(:,:),1),1);
    elseif strcmp(sessionInfo.subject, 'MrJones')
        %jones
        low_current = [0.5 NaN];
        high_current = [1.0 1.5];
%         tmp_trialinfo(:,5) = ones(size(ALL_TFRbs{1}.trialinfo(:,:),1),1);
    end
    
    %Hi vs low current
%     current_levels = [0.5 1.0 1.5 2.0 2.5];
    current_trials_low = find(data.trialinfo(:,2) == low_current(1) | data.trialinfo(:,2) == low_current(2));
    current_trials_high = find(data.trialinfo(:,2) == high_current(1) | data.trialinfo(:,2) == high_current(2));
    
    figure('Position', [25 10 2500 980], 'Name', ['Unit DBS Analysis, Session = ' int2str(s)])
    %PFC
    subplot(3,4,1);
    imagesc( squeeze(nanmean(spike_rates{1}(current_trials_low,:,:),1))); 
    colormap('Jet')
    caxis([0 5])
    colorbar
    title(['PFC - low current, All chans (abs)'])
    
    subplot(3,4,5);
    imagesc( squeeze(nanmean(spike_rates{1}(current_trials_high,:,:),1))); 
    colormap('Jet')
    caxis([0 5])
    colorbar
    title(['PFC - high current, All chans (abs)'])    
    
    subplot(3,4,9);
    plot(time_axis, squeeze(nanmean(nanmean(spike_rates_bs{1}(current_trials_high,:,:),1),2)), 'r'); hold on
    plot(time_axis, squeeze(nanmean(nanmean(spike_rates_bs{1}(current_trials_low,:,:),1),2)), 'b')
    title(['low vs. high current, All chans (bs)'])

    %FEF
    subplot(3,4,2);
    imagesc( squeeze(nanmean(spike_rates{2}(current_trials_low,:,:),1))); 
    colormap('Jet')
    caxis([0 5])
    colorbar
    title(['FEF - low current, All chans (abs)'])
    
    subplot(3,4,6);
    imagesc( squeeze(nanmean(spike_rates{2}(current_trials_high,:,:),1))); 
    colormap('Jet')
    caxis([0 5])
    colorbar
    title(['FEF - high current, All chans (abs)'])    
    
    subplot(3,4,10);
    plot(time_axis, squeeze(nanmean(nanmean(spike_rates_bs{2}(current_trials_high,:,:),1),2)), 'r'); hold on
    plot(time_axis, squeeze(nanmean(nanmean(spike_rates_bs{2}(current_trials_low,:,:),1),2)), 'b')
    title(['low vs. high current, All chans (bs)'])

    %PPC
    subplot(3,4,3);
    imagesc( squeeze(nanmean(spike_rates{3}(current_trials_low,:,:),1))); 
    colormap('Jet')
    caxis([0 5])
    colorbar
    title(['PPC - low current, All chans (abs)'])
    
    subplot(3,4,7);
    imagesc( squeeze(nanmean(spike_rates{3}(current_trials_high,:,:),1))); 
    colormap('Jet')
    caxis([0 5])
    colorbar
    title(['PPC - high current, All chans (abs)'])    
    
    subplot(3,4,11);
    plot(time_axis, squeeze(nanmean(nanmean(spike_rates_bs{3}(current_trials_high,:,:),1),2)), 'r'); hold on
    plot(time_axis, squeeze(nanmean(nanmean(spike_rates_bs{3}(current_trials_low,:,:),1),2)), 'b')
    title(['low vs. high current, All chans (bs)'])

    %STG
    subplot(3,4,4);
    imagesc( squeeze(nanmean(spike_rates{4}(current_trials_low,:,:),1))); 
    colormap('Jet')
    caxis([0 5])
    colorbar
    title(['STG - low current, All chans (abs)'])
    
    subplot(3,4,8);
    imagesc( squeeze(nanmean(spike_rates{4}(current_trials_high,:,:),1))); 
    colormap('Jet')
    caxis([0 5])
    colorbar
    title(['STG - high current, All chans (abs)'])    
    
    subplot(3,4,12);
    plot(time_axis, squeeze(nanmean(nanmean(spike_rates_bs{4}(current_trials_high,:,:),1),2)), 'r'); hold on
    plot(time_axis, squeeze(nanmean(nanmean(spike_rates_bs{4}(current_trials_low,:,:),1),2)), 'b')
    title(['low vs. high current, All chans (bs)'])
    
    
    cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')
    save([sessionInfo.session '_DBS_sortedSUA'], 'spike_rates_bs', 'spike_rates', 'current_trials_high', 'current_trials_low', ...
        'sessionInfo','trialInfo', 'kept_vs_all', ...
        'time_axis', 'newTrialInfo', 'standardTrialIndx',  '-v7.3')
    
    clear spike_rates_bs spike_rates current_trials_high current_trials_low sessionInfo trialInfo newTrialInfo standardTrialIndx
end



return


%%
filestoload = ls('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\*_DBS_sortedSUA*.mat');


current_levels = [0.5 1.0 1.5 2.0 2.5];
monkey_index = {};
spike_rates_bs_sessions = {};
spike_rates_sessions = {};
trialInfo_sessions = [];
for s = 1:size(filestoload,1)
    cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')
    load(deblank(filestoload(s,:)))
    
    if strcmp(sessionInfo.subject, 'Mary')
        %mary
        low_current = [0.5 1.0];
        high_current = [1.5 2.5];
    elseif strcmp(sessionInfo.subject, 'MrJones')
        %jones
        low_current = [0.5 NaN];
        high_current = [1.0 1.5];
    end
        
    %low current
    current_trials_low = find((newTrialInfo(:,2) == low_current(1) | newTrialInfo(:,2) == low_current(2)) );
    current_trials_high = find((newTrialInfo(:,2) == high_current(1) | newTrialInfo(:,2) == high_current(2)) );    
    
    if s == 1
        %low
        spike_rates_bs_sessions{1,1} = squeeze(mean(spike_rates{1}(current_trials_low,:,:),1));
        spike_rates_bs_sessions{2,1} = squeeze(mean(spike_rates{2}(current_trials_low,:,:),1));
        spike_rates_bs_sessions{3,1} = squeeze(mean(spike_rates{3}(current_trials_low,:,:),1));
%         spike_rates_bs_sessions{4,1} =
%         squeeze(mean(spike_rates_bs{4}(current_trials_low,:,:),1)); %need
%         to do this weirdness for first session because only 1 spikeing
%         neuron in STG
        spike_rates_bs_sessions{4,1} = squeeze(mean(spike_rates{4}(current_trials_low,:,:),1))';
        %high
        spike_rates_bs_sessions{1,2} = squeeze(mean(spike_rates{1}(current_trials_high,:,:),1));
        spike_rates_bs_sessions{2,2} = squeeze(mean(spike_rates{2}(current_trials_high,:,:),1));
        spike_rates_bs_sessions{3,2} = squeeze(mean(spike_rates{3}(current_trials_high,:,:),1));
%         spike_rates_bs_sessions{4,2} = squeeze(mean(spike_rates_bs{4}(current_trials_high,:,:),1));        
        spike_rates_bs_sessions{4,2} = squeeze(mean(spike_rates{4}(current_trials_high,:,:),1))';  
        
        %save monkey info
        if strcmp(sessionInfo.subject, 'Mary')
            monkey_index{1} =  2.*ones(size(spike_rates{1},2),1);
            monkey_index{2} =  2.*ones(size(spike_rates{2},2),1);
            monkey_index{3} =  2.*ones(size(spike_rates{3},2),1);
            monkey_index{4} =  2.*ones(size(spike_rates{4},2),1);
        elseif strcmp(sessionInfo.subject, 'MrJones')
            monkey_index{1} =  ones(size(spike_rates{1},2),1);
            monkey_index{2} =  ones(size(spike_rates{2},2),1);
            monkey_index{3} =  ones(size(spike_rates{3},2),1);
            monkey_index{4} =  ones(size(spike_rates{4},2),1);
        end
        trialInfo_sessions = newTrialInfo;
    else

        %low
        spike_rates_bs_sessions{1,1} = cat(1, spike_rates_bs_sessions{1,1}, squeeze(mean(spike_rates{1}(current_trials_low,:,:),1)));
        spike_rates_bs_sessions{2,1} = cat(1, spike_rates_bs_sessions{2,1}, squeeze(mean(spike_rates{2}(current_trials_low,:,:),1)));
        spike_rates_bs_sessions{3,1} = cat(1, spike_rates_bs_sessions{3,1}, squeeze(mean(spike_rates{3}(current_trials_low,:,:),1)));
        if ~isempty(spike_rates{4})
            spike_rates_bs_sessions{4,1} = cat(1, spike_rates_bs_sessions{4,1}, squeeze(mean(spike_rates{4}(current_trials_low,:,:),1)));
        end

        %high
        spike_rates_bs_sessions{1,2} = cat(1, spike_rates_bs_sessions{1,2}, squeeze(mean(spike_rates{1}(current_trials_high,:,:),1)));
        spike_rates_bs_sessions{2,2} = cat(1, spike_rates_bs_sessions{2,2}, squeeze(mean(spike_rates{2}(current_trials_high,:,:),1)));
        spike_rates_bs_sessions{3,2} = cat(1, spike_rates_bs_sessions{3,2}, squeeze(mean(spike_rates{3}(current_trials_high,:,:),1)));
        if ~isempty(spike_rates{4})
            spike_rates_bs_sessions{4,2} = cat(1, spike_rates_bs_sessions{4,2}, squeeze(mean(spike_rates{4}(current_trials_high,:,:),1)));
        end
%         spike_rates_sessions{1} = cat(1, spike_rates_sessions{1}, spike_rates{1});
%         spike_rates_sessions{2} = cat(1, spike_rates_sessions{2}, spike_rates{2});
%         spike_rates_sessions{3} = cat(1, spike_rates_sessions{3}, spike_rates{3});
%         spike_rates_sessions{4} = cat(1, spike_rates_sessions{4}, spike_rates{4});

        %save monkey info
        if strcmp(sessionInfo.subject, 'Mary')
            if ~isempty(monkey_index{1}), monkey_index{1} =  [monkey_index{1};  2.*ones(size(spike_rates{1},2),1)]; end
            if ~isempty(monkey_index{2}), monkey_index{2} =  [monkey_index{2};  2.*ones(size(spike_rates{2},2),1)]; end
            if ~isempty(monkey_index{3}), monkey_index{3} =  [monkey_index{3};  2.*ones(size(spike_rates{3},2),1)]; end
            if ~isempty(monkey_index{4}), monkey_index{4} =  [monkey_index{4};  2.*ones(size(spike_rates{4},2),1)]; end
        elseif strcmp(sessionInfo.subject, 'MrJones')
            if ~isempty(monkey_index{1}), monkey_index{1} =  [monkey_index{1};  ones(size(spike_rates{1},2),1)]; end
            if ~isempty(monkey_index{1}), monkey_index{2} =  [monkey_index{2};  ones(size(spike_rates{2},2),1)]; end
            if ~isempty(monkey_index{1}), monkey_index{3} =  [monkey_index{3};  ones(size(spike_rates{3},2),1)]; end
            if ~isempty(monkey_index{1}), monkey_index{4} =  [monkey_index{4};  ones(size(spike_rates{4},2),1)]; end
        end
        trialInfo_sessions = [trialInfo_sessions; newTrialInfo];
    end
    

    areas = {'PFC', '8A', 'PPC', 'STG'};

   
end


nunits1 = size(spike_rates_bs_sessions{1,1},1);
nunits2 = size(spike_rates_bs_sessions{1,2},1);
sem1 = nanstd(spike_rates_bs_sessions{1,1})./ sqrt(nunits1);
sem2 = nanstd(spike_rates_bs_sessions{1,2})./ sqrt(nunits2);
figure; shadedErrorBar(time_axis, nanmean(spike_rates_bs_sessions{1,1}), 2.*sem1, 'b');
hold on; shadedErrorBar(time_axis, nanmean(spike_rates_bs_sessions{1,2}), 2.*sem2,  'r');
yline(0, 'k--')
xlabel('Time from electrical stim onset (s)')
ylabel('Spike rate change from baseline (spikes/sec)')
xlim([-30 118])
title(['PFC N = ' int2str(nunits2)])

nunits1 = size(spike_rates_bs_sessions{2,1},1);
nunits2 = size(spike_rates_bs_sessions{2,2},1);
sem1 = nanstd(spike_rates_bs_sessions{2,1})./ sqrt(nunits1);
sem2 = nanstd(spike_rates_bs_sessions{2,2})./ sqrt(nunits2);
figure; shadedErrorBar(time_axis, nanmean(spike_rates_bs_sessions{2,1}), 2.*sem1, 'b');
hold on; shadedErrorBar(time_axis, nanmean(spike_rates_bs_sessions{2,2}), 2.*sem2, 'r');
yline(0, 'k--')
xlabel('Time from electrical stim onset (s)')
ylabel('Spike rate change from baseline (spikes/sec)')
xlim([-30 118])
title(['FEF N= ' int2str(nunits2)])

nunits1 = size(spike_rates_bs_sessions{3,1},1);
nunits2 = size(spike_rates_bs_sessions{3,2},1);
sem1 = nanstd(spike_rates_bs_sessions{3,1})./ sqrt(nunits1);
sem2 = nanstd(spike_rates_bs_sessions{3,2})./ sqrt(nunits2);
figure; shadedErrorBar(time_axis, nanmean(spike_rates_bs_sessions{3,1}), 2.*sem1, 'b');
hold on; shadedErrorBar(time_axis, nanmean(spike_rates_bs_sessions{3,2}), 2.*sem2, 'r');
yline(0, 'k--')
xlabel('Time from electrical stim onset (s)')
ylabel('Spike rate change from baseline (spikes/sec)')
xlim([-30 118])
title(['PPC N= ' int2str(nunits2)])

nunits1 = size(spike_rates_bs_sessions{4,1},1);
nunits2 = size(spike_rates_bs_sessions{4,2},1);
sem1 = nanstd(spike_rates_bs_sessions{4,1})./ sqrt(nunits1);
sem2 = nanstd(spike_rates_bs_sessions{4,2})./ sqrt(nunits2);
figure; shadedErrorBar(time_axis, nanmean(spike_rates_bs_sessions{4,1}), 2.*sem1, 'b');
hold on; shadedErrorBar(time_axis, nanmean(spike_rates_bs_sessions{4,2}), 2.*sem2, 'r');
yline(0, 'k--')
xlabel('Time from electrical stim onset (s)')
ylabel('Spike rate change from baseline (spikes/sec)')
xlim([-30 118])
title(['STG N= ' int2str(nunits2)])


    %create time axis
    xmarkers = 1:5:length(time_axis);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis(xmarkers(indx)));
        indx=indx+1;
    end
    
figure; imagesc(spike_rates_bs_sessions{4,2}./max(spike_rates_bs_sessions{4,2})); caxis([0 1])

figure; 
subplot(1,2,1);
imagesc(spike_rates_bs_sessions{1,1}); caxis([0 5]); 
hold on; plot(25.*ones(length(find(monkey_index{1}==1)),1), find(monkey_index{1}==1), 'w')
set(gca, 'xtick', xmarkers)
set(gca, 'xticklabel', xlabels)
title('PFC Low current')

subplot(1,2,2);
imagesc(spike_rates_bs_sessions{1,2}); caxis([0 5])
hold on; plot(25.*ones(length(find(monkey_index{1}==1)),1), find(monkey_index{1}==1), 'w')
set(gca, 'xtick', xmarkers)
set(gca, 'xticklabel', xlabels)
title('PFC High current')

figure; 
subplot(1,2,1);
imagesc(spike_rates_bs_sessions{2,1}); caxis([0 5]); 
hold on; plot(25.*ones(length(find(monkey_index{2}==1)),1), find(monkey_index{2}==1), 'w')
set(gca, 'xtick', xmarkers)
set(gca, 'xticklabel', xlabels)
title('8A Low current')

subplot(1,2,2);
imagesc(spike_rates_bs_sessions{2,2}); caxis([0 5])
hold on; plot(25.*ones(length(find(monkey_index{2}==1)),1), find(monkey_index{2}==1), 'w')
set(gca, 'xtick', xmarkers)
set(gca, 'xticklabel', xlabels)
title('8A High current')


figure; 
subplot(1,2,1);
imagesc(spike_rates_bs_sessions{4,1}); caxis([0 5]); 
hold on; plot(25.*ones(length(find(monkey_index{4}==1)),1), find(monkey_index{4}==1), 'w')
set(gca, 'xtick', xmarkers)
set(gca, 'xticklabel', xlabels)
title('STG Low current')

subplot(1,2,2);
imagesc(spike_rates_bs_sessions{4,2}); caxis([0 5])
hold on; plot(25.*ones(length(find(monkey_index{4}==1)),1), find(monkey_index{4}==1), 'w')
set(gca, 'xtick', xmarkers)
set(gca, 'xticklabel', xlabels)
title('STG High current')

figure; 
subplot(1,2,1);
imagesc(spike_rates_bs_sessions{3,1}); caxis([0 5]); 
hold on; plot(25.*ones(length(find(monkey_index{3}==1)),1), find(monkey_index{3}==1), 'w')
set(gca, 'xtick', xmarkers)
set(gca, 'xticklabel', xlabels)
title('PPC Low current')

subplot(1,2,2);
imagesc(spike_rates_bs_sessions{3,2}); caxis([0 5])
hold on; plot(25.*ones(length(find(monkey_index{3}==1)),1), find(monkey_index{3}==1), 'w')
set(gca, 'xtick', xmarkers)
set(gca, 'xticklabel', xlabels)
title('PPC High current')

    
all_spikes_high = cat(1, spike_rates_bs_sessions{1,2}, spike_rates_bs_sessions{2,2}, spike_rates_bs_sessions{3,2}, spike_rates_bs_sessions{4,2});

figure; imagesc(all_spikes_high); caxis([0 4])
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
colorbar
xlabel('Time relative to DBS onset (sec)');
ylabel('Unit index')
title('Firing Rates during high current stimulation')

all_spikes_low = cat(1, spike_rates_bs_sessions{1,1}, spike_rates_bs_sessions{2,1}, spike_rates_bs_sessions{3,1}, spike_rates_bs_sessions{4,1});

figure; imagesc(all_spikes_low); caxis([0 4])
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
colorbar
xlabel('Time relative to DBS onset (sec)');
ylabel('Unit index')
title('Firing Rates during low current stimulation')




nunits1 = size(all_spikes_high,1);
nunits2 = size(all_spikes_low,1);
sem1 = nanstd(all_spikes_high)./ sqrt(nunits1);
sem2 = nanstd(all_spikes_low)./ sqrt(nunits2);
figure; 
hold on; shadedErrorBar(time_axis, nanmean(all_spikes_high), 2.*sem1, 'r');
shadedErrorBar(time_axis, nanmean(all_spikes_low), 2.*sem2, 'b');
xline(0, 'k--')
xline(nanmean(trialInfo_sessions(:,4)), 'k--')
xlabel('Time from electrical stim onset (s)')
ylabel('Firing rate (spikes/sec)')
xlim([-30 120])
set(gca, 'xtick', [-30:30:120])
title(['All Units, N = ' int2str(size(all_spikes_high,1)) ', Blue = Low Current, Red = High Current'])




%%
%With SEM and ttest
num_neurons1 = size(spike_rates_bs_sessions{1,1},1);
num_neurons2 = size(spike_rates_bs_sessions{2,1},1);
num_neurons3 = size(spike_rates_bs_sessions{3,1},1);
num_neurons4 = size(spike_rates_bs_sessions{4,1},1);

%pre-stim
t1 = nearest(time_axis, -30);
t2 = nearest(time_axis, -1);
% spike_rates_bs_sessions{1,1} = spike_rates_bs_sessions{1,1} -  nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2),2);
% spike_rates_bs_sessions{2,1} = spike_rates_bs_sessions{2,1} -  nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2),2);
% spike_rates_bs_sessions{3,1} = spike_rates_bs_sessions{3,1} -  nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2),2);
% spike_rates_bs_sessions{4,1} = spike_rates_bs_sessions{4,1} -  nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2),2);
% spike_rates_bs_sessions{1,2} = spike_rates_bs_sessions{1,2} -  nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2);
% spike_rates_bs_sessions{2,2} = spike_rates_bs_sessions{2,2} -  nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2);
% spike_rates_bs_sessions{3,2} = spike_rates_bs_sessions{3,2} -  nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2);
% spike_rates_bs_sessions{4,2} = spike_rates_bs_sessions{4,2} -  nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2);
%Bar plots
spikes_bar = zeros(4,5,2); %Dim is area x time period ( pre-stim, stim, post1/2/3) x current (low,high)
spikes_bar_sem = zeros(4,5,2);
spikes_stats = zeros(4,5,2);
%pre-stim
bs1 = nearest(time_axis, -30);
bs2 = nearest(time_axis, -1);
%Low Current
spikes_bar(1,1,1) = nanmean(nanmean(spike_rates_bs_sessions{1,1}(:,bs1:bs2)));
spikes_bar(2,1,1) = nanmean(nanmean(spike_rates_bs_sessions{2,1}(:,bs1:bs2)));
spikes_bar(3,1,1) = nanmean(nanmean(spike_rates_bs_sessions{3,1}(:,bs1:bs2)));
spikes_bar(4,1,1) = nanmean(nanmean(spike_rates_bs_sessions{4,1}(:,bs1:bs2)));
%High Current
spikes_bar(1,1,2) = nanmean(nanmean(spike_rates_bs_sessions{1,2}(:,bs1:bs2)));
spikes_bar(2,1,2) = nanmean(nanmean(spike_rates_bs_sessions{2,2}(:,bs1:bs2)));
spikes_bar(3,1,2) = nanmean(nanmean(spike_rates_bs_sessions{3,2}(:,bs1:bs2)));
spikes_bar(4,1,2) = nanmean(nanmean(spike_rates_bs_sessions{4,2}(:,bs1:bs2)));
%Low current SEM
spikes_bar_sem(1,1,1) = nanstd(nanmean(spike_rates_bs_sessions{1,1}(:,bs1:bs2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,1,1) = nanstd(nanmean(spike_rates_bs_sessions{2,1}(:,bs1:bs2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,1,1) = nanstd(nanmean(spike_rates_bs_sessions{3,1}(:,bs1:bs2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,1,1) = nanstd(nanmean(spike_rates_bs_sessions{4,1}(:,bs1:bs2),2))./sqrt(num_neurons4);
%High current SEM
spikes_bar_sem(1,1,2) = nanstd(nanmean(spike_rates_bs_sessions{1,2}(:,bs1:bs2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,1,2) = nanstd(nanmean(spike_rates_bs_sessions{2,2}(:,bs1:bs2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,1,2) = nanstd(nanmean(spike_rates_bs_sessions{3,2}(:,bs1:bs2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,1,2) = nanstd(nanmean(spike_rates_bs_sessions{4,2}(:,bs1:bs2),2))./sqrt(num_neurons4);
%stim
t1 = nearest(time_axis, 1);
t2 = nearest(time_axis, 28);
%Low Current
spikes_bar(1,2,1) = nanmean(nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2)));
spikes_bar(2,2,1) = nanmean(nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2)));
spikes_bar(3,2,1) = nanmean(nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2)));
spikes_bar(4,2,1) = nanmean(nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2)));
%High Current
spikes_bar(1,2,2) = nanmean(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2)));
spikes_bar(2,2,2) = nanmean(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2)));
spikes_bar(3,2,2) = nanmean(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2)));
spikes_bar(4,2,2) = nanmean(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2)));
%Low current SEM
spikes_bar_sem(1,2,1) = nanstd(nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,2,1) = nanstd(nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,2,1) = nanstd(nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,2,1) = nanstd(nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2),2))./sqrt(num_neurons4);
%High current SEM
spikes_bar_sem(1,2,2) = nanstd(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,2,2) = nanstd(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,2,2) = nanstd(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,2,2) = nanstd(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2))./sqrt(num_neurons4);
%stat test
for a = 1:4
  for c = 1:2
    tmp1 = log(nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2));
    tmp2 = log(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2));
    rmi = find(tmp1<-10E10 | tmp2 <-10E10);
    tmp1(rmi) = []; tmp2(rmi) = [];
    [h p] = ttest2(tmp1, tmp2);
    spikes_stats(a,2,c) = p;
  end
end

%post-stim1
t1 = nearest(time_axis, 31);
t2 = nearest(time_axis, 60);
%Low Current
spikes_bar(1,3,1) = nanmean(nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2)));
spikes_bar(2,3,1) = nanmean(nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2)));
spikes_bar(3,3,1) = nanmean(nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2)));
spikes_bar(4,3,1) = nanmean(nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2)));
%High Current
spikes_bar(1,3,2) = nanmean(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2)));
spikes_bar(2,3,2) = nanmean(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2)));
spikes_bar(3,3,2) = nanmean(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2)));
spikes_bar(4,3,2) = nanmean(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2)));
%Low current SEM
spikes_bar_sem(1,3,1) = nanstd(nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,3,1) = nanstd(nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,3,1) = nanstd(nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,3,1) = nanstd(nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2),2))./sqrt(num_neurons4);
%High current SEM
spikes_bar_sem(1,3,2) = nanstd(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,3,2) = nanstd(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,3,2) = nanstd(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,3,2) = nanstd(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2))./sqrt(num_neurons4);
%stat test
for a = 1:4
  for c = 1:2
    tmp1 = log(nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2));
    tmp2 = log(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2));
    rmi = find(tmp1<-10E10 | tmp2 <-10E10);
    tmp1(rmi) = []; tmp2(rmi) = [];
    [h p] = ttest2(tmp1, tmp2);
    spikes_stats(a,3,c) = p;
  end
end
%post-stim2
t1 = nearest(time_axis, 61);
t2 = nearest(time_axis, 90);
%Low Current
spikes_bar(1,4,1) = nanmean(nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2)));
spikes_bar(2,4,1) = nanmean(nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2)));
spikes_bar(3,4,1) = nanmean(nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2)));
spikes_bar(4,4,1) = nanmean(nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2)));
%High Current
spikes_bar(1,4,2) = nanmean(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2)));
spikes_bar(2,4,2) = nanmean(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2)));
spikes_bar(3,4,2) = nanmean(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2)));
spikes_bar(4,4,2) = nanmean(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2)));
%Low current SEM
spikes_bar_sem(1,4,1) = nanstd(nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,4,1) = nanstd(nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,4,1) = nanstd(nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,4,1) = nanstd(nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2),2))./sqrt(num_neurons4);
%High current SEM
spikes_bar_sem(1,4,2) = nanstd(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,4,2) = nanstd(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,4,2) = nanstd(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,4,2) = nanstd(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2))./sqrt(num_neurons4);
%stat test
for a = 1:4
  for c = 1:2
    tmp1 = log(nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2));
    tmp2 = log(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2));
    rmi = find(tmp1<-10E10 | tmp2 <-10E10);
    tmp1(rmi) = []; tmp2(rmi) = [];
    [h p] = ttest2(tmp1, tmp2);
    spikes_stats(a,4,c) = p;
  end
end

%post-stim3
t1 = nearest(time_axis, 91);
t2 = nearest(time_axis, 120);
%Low Current
spikes_bar(1,5,1) = nanmean(nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2),2));
spikes_bar(2,5,1) = nanmean(nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2),2));
spikes_bar(3,5,1) = nanmean(nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2),2));
spikes_bar(4,5,1) = nanmean(nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2),2));
%High Current
spikes_bar(1,5,2) = nanmean(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2));
spikes_bar(2,5,2) = nanmean(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2));
spikes_bar(3,5,2) = nanmean(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2));
spikes_bar(4,5,2) = nanmean(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2));
%Low current SEM
spikes_bar_sem(1,5,1) = nanstd(nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,5,1) = nanstd(nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,5,1) = nanstd(nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,5,1) = nanstd(nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2),2))./sqrt(num_neurons4);
%High current SEM
spikes_bar_sem(1,5,2) = nanstd(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2))./sqrt(num_neurons1);
spikes_bar_sem(2,5,2) = nanstd(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2))./sqrt(num_neurons2);
spikes_bar_sem(3,5,2) = nanstd(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2))./sqrt(num_neurons3);
spikes_bar_sem(4,5,2) = nanstd(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2))./sqrt(num_neurons4);
%stat test
for a = 1:4
  for c = 1:2
    tmp1 = log(nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2));
    tmp2 = log(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2));
    rmi = find(tmp1<-10E10 | tmp2 <-10E10);
    tmp1(rmi) = []; tmp2(rmi) = [];
    [h p] = ttest2(tmp1, tmp2);
    spikes_stats(a,5,c) = p;
  end
end
figure;
subplot(2,1,1);
bar(squeeze(spikes_bar(:,:,2))); hold on
errorbar([1:4]-0.31, squeeze(spikes_bar(:,1,2)), 2.*squeeze(spikes_bar_sem(:,1,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]-0.16, squeeze(spikes_bar(:,2,2)), 2.*squeeze(spikes_bar_sem(:,2,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4], squeeze(spikes_bar(:,3,2)), 2.*squeeze(spikes_bar_sem(:,3,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.16, squeeze(spikes_bar(:,4,2)), 2.*squeeze(spikes_bar_sem(:,4,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.31, squeeze(spikes_bar(:,5,2)), 2.*squeeze(spikes_bar_sem(:,5,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
xpos = [NaN -0.16 0 0.16 0.31];
for a = 1:4
   for tp = 2:5
       if spikes_stats(a,tp,2)<0.01
           plot(a+xpos(tp), squeeze(spikes_bar(a,tp,2))*1.3, 'r*')
       end
   end
end
%Firing rates during natural ROC
plot([0.5 1.5], [3.4555 3.4555], 'b--'); %PFC
plot([1.5 2.5], [3.1227 3.1227], 'r--'); %8A
plot([2.5 3.5], [3.5027 3.5027], 'g--'); %PPC
plot([3.5 4.5], [5.0097 5.0097], 'm--'); %STG
ylim([0.5 7])
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'})
ylabel('Firing Rate (spikes/sec)')
title('High Current')
subplot(2,1,2);
bar(squeeze(spikes_bar(:,:,1))); hold on;
errorbar([1:4]-0.31, squeeze(spikes_bar(:,1,1)), 2.*squeeze(spikes_bar_sem(:,1,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]-0.16, squeeze(spikes_bar(:,2,1)), 2.*squeeze(spikes_bar_sem(:,2,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.00, squeeze(spikes_bar(:,3,1)), 2.*squeeze(spikes_bar_sem(:,3,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.16, squeeze(spikes_bar(:,4,1)), 2.*squeeze(spikes_bar_sem(:,4,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.31, squeeze(spikes_bar(:,5,1)), 2.*squeeze(spikes_bar_sem(:,5,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
for a = 1:4
   for tp = 2:5
       if spikes_stats(a,tp,1)<0.01
           plot(a+xpos(tp), squeeze(spikes_bar(a,tp,1))*1.3, 'r*')
       end
   end
end
%Firing rates during natural ROC
plot([0.5 1.5], [3.4555 3.4555], 'b--'); %PFC
plot([1.5 2.5], [3.1227 3.1227], 'r--'); %8A
plot([2.5 3.5], [3.5027 3.5027], 'g--'); %PPC
plot([3.5 4.5], [5.0097 5.0097], 'm--'); %STG
ylim([0.5 7])
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'})
ylabel('Firing Rate (spikes/sec)')
title('Low Current')


t1 = nearest(time_axis, 1);
t2 = nearest(time_axis, 28);

[h p] = signtest(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2)-5.0097)
[h p] = signtest(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2)-3.5)
[h p] = signtest(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2)-3.12)
[h p] = signtest(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2)-3.45)


[h p] = ttest(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2)-5.0097)


[h p] = ttest(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2)-3.5)
[h p] = ttest(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2)-3.12)
[h p] = ttest(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2)-3.45)




%%
addpath('C:\Codes\')

%With nonparametric CI and randomization test
num_neurons1 = size(spike_rates_bs_sessions{1,1},1);
num_neurons2 = size(spike_rates_bs_sessions{2,1},1);
num_neurons3 = size(spike_rates_bs_sessions{3,1},1);
num_neurons4 = size(spike_rates_bs_sessions{4,1},1);

%pre-stim
t1 = nearest(time_axis, -30);
t2 = nearest(time_axis, -1);
% spike_rates_bs_sessions{1,1} = spike_rates_bs_sessions{1,1} -  nanmean(spike_rates_bs_sessions{1,1}(:,t1:t2),2);
% spike_rates_bs_sessions{2,1} = spike_rates_bs_sessions{2,1} -  nanmean(spike_rates_bs_sessions{2,1}(:,t1:t2),2);
% spike_rates_bs_sessions{3,1} = spike_rates_bs_sessions{3,1} -  nanmean(spike_rates_bs_sessions{3,1}(:,t1:t2),2);
% spike_rates_bs_sessions{4,1} = spike_rates_bs_sessions{4,1} -  nanmean(spike_rates_bs_sessions{4,1}(:,t1:t2),2);
% spike_rates_bs_sessions{1,2} = spike_rates_bs_sessions{1,2} -  nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2);
% spike_rates_bs_sessions{2,2} = spike_rates_bs_sessions{2,2} -  nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2);
% spike_rates_bs_sessions{3,2} = spike_rates_bs_sessions{3,2} -  nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2);
% spike_rates_bs_sessions{4,2} = spike_rates_bs_sessions{4,2} -  nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2);
%Bar plots
spikes_bar = zeros(4,5,2); %Dim is area x time period ( pre-stim, stim, post1/2/3) x current (low,high)
% spikes_CI_H = zeros(4,5,2);
% spikes_CI_L = zeros(4,5,2);
spikes_stats = zeros(4,5,2);
%pre-stim
bs1 = nearest(time_axis, -30);
bs2 = nearest(time_axis, -1);
for a =1:4
  for c = 1:2 %low/high current
    spikes_bar(a,1,c) = nanmean(nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2)));
    spikes_bar(a,1,c)
    [confInts,statResmp,confIdxs]  = bootstrapConfIntervals(nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2), @nanmean, [],0.99,10000);
    spikes_CI_H(a,1,c) = confInts(1);
    spikes_CI_L(a,1,c) = confInts(2);
  end
%   a
end
%stim
t1 = nearest(time_axis, 1);
t2 = nearest(time_axis, 28);
for a =1:4
  for c = 1:2 %low/high current
    spikes_bar(a,2,c) = nanmean(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2)));
    [confInts,statResmp,confIdxs]  = bootstrapConfIntervals(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2), @nanmean, [],0.99,10000);
    spikes_CI_H(a,2,c) = confInts(1);
    spikes_CI_L(a,2,c) = confInts(2);
  end
  a
end

%stat test
for a = 1:4
  for c = 1:2
    tmp1 = (nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2));
    tmp2 = (nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2));
%     rmi = find(tmp1<-10E10 | tmp2 <-10E10);
%     tmp1(rmi) = []; tmp2(rmi) = [];
    [p h] = nonparam_mean_diff(tmp1, tmp2);
    spikes_stats(a,2,c) = p;
    title(['stim period, area: ' int2str(a) ' current: ' int2str(c)])
  end
end

%post-stim1
t1 = nearest(time_axis, 31);
t2 = nearest(time_axis, 60);
for a =1:4
  for c = 1:2 %low/high current
    spikes_bar(a,3,c) = nanmean(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2)));
    [confInts,statResmp,confIdxs]  = bootstrapConfIntervals(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2), @nanmean, [],0.99,10000);
    spikes_CI_H(a,3,c) = confInts(1);
    spikes_CI_L(a,3,c) = confInts(2);
  end
  a
end
%stat test
for a = 1:4
  for c = 1:2
    tmp1 = (nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2));
    tmp2 = (nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2));
%     rmi = find(tmp1<-10E10 | tmp2 <-10E10);
%     tmp1(rmi) = []; tmp2(rmi) = [];
    [p h] = nonparam_mean_diff(tmp1, tmp2);
    spikes_stats(a,3,c) = p;
    title(['post1 period, area: ' int2str(a) ' current: ' int2str(c)])    
  end
end
%post-stim2
t1 = nearest(time_axis, 61);
t2 = nearest(time_axis, 90);
for a =1:4
  for c = 1:2 %low/high current
    spikes_bar(a,4,c) = nanmean(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2)));
    [confInts,statResmp,confIdxs]  = bootstrapConfIntervals(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2), @nanmean, [],0.99,10000);
    spikes_CI_H(a,4,c) = confInts(1);
    spikes_CI_L(a,4,c) = confInts(2);
  end
  a
end
%stat test
for a = 1:4
  for c = 1:2
    tmp1 = (nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2));
    tmp2 = (nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2));
%     rmi = find(tmp1<-10E10 | tmp2 <-10E10);
%     tmp1(rmi) = []; tmp2(rmi) = [];
    [p h] = nonparam_mean_diff(tmp1, tmp2);
    spikes_stats(a,4,c) = p;
    title(['post2 period, area: ' int2str(a) ' current: ' int2str(c)])        
  end
end

%post-stim3
t1 = nearest(time_axis, 91);
t2 = nearest(time_axis, 120);
for a =1:4
  for c = 1:2 %low/high current
    spikes_bar(a,5,c) = nanmean(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2)));
    [confInts,statResmp,confIdxs]  = bootstrapConfIntervals(nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2), @nanmean, [],0.99,10000);
    spikes_CI_H(a,5,c) = confInts(1);
    spikes_CI_L(a,5,c) = confInts(2);
  end
  a
end
%stat test
for a = 1:4
  for c = 1:2
    tmp1 = (nanmean(spike_rates_bs_sessions{a,c}(:,bs1:bs2),2));
    tmp2 = (nanmean(spike_rates_bs_sessions{a,c}(:,t1:t2),2));
%     rmi = find(tmp1<-10E10 | tmp2 <-10E10);
%     tmp1(rmi) = []; tmp2(rmi) = [];
    [p h] = nonparam_mean_diff(tmp1, tmp2);
    if length(p)>1
        p=p(1);
    end
    spikes_stats(a,5,c) = p;
    title(['post3 period, area: ' int2str(a) ' current: ' int2str(c)])    
    
  end
end

close all

return 

%%
figure;
subplot(2,1,1);
bar(squeeze(spikes_bar(:,:,2))); hold on
%pre-stim
NEG = squeeze(spikes_bar(:,1,2)-spikes_CI_L(:,1,2)); POS = squeeze(spikes_CI_H(:,1,2)-squeeze(spikes_bar(:,1,2)));
errorbar([1:4]-0.31, squeeze(spikes_bar(:,1,2)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)
%stim
NEG = squeeze(spikes_bar(:,2,2)-spikes_CI_L(:,2,2)); POS = squeeze(spikes_CI_H(:,1,2)-squeeze(spikes_bar(:,1,2)));
errorbar([1:4]-0.16, squeeze(spikes_bar(:,2,2)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)
%post1
NEG = squeeze(spikes_bar(:,3,2)-spikes_CI_L(:,3,2)); POS = squeeze(spikes_CI_H(:,1,2)-squeeze(spikes_bar(:,1,2)));
errorbar([1:4], squeeze(spikes_bar(:,3,2)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)
%post2
NEG = squeeze(spikes_bar(:,4,2)-spikes_CI_L(:,4,2)); POS = squeeze(spikes_CI_H(:,1,2)-squeeze(spikes_bar(:,1,2)));
errorbar([1:4]+0.16, squeeze(spikes_bar(:,4,2)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)
%post3
NEG = squeeze(spikes_bar(:,5,2)-spikes_CI_L(:,5,2)); POS = squeeze(spikes_CI_H(:,1,2)-squeeze(spikes_bar(:,1,2)));
errorbar([1:4]+0.31, squeeze(spikes_bar(:,5,2)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)

xpos = [NaN -0.16 0 0.16 0.31];
for a = 1:4
   for tp = 2:5
       if spikes_stats(a,tp,2)<0.01
           plot(a+xpos(tp), squeeze(spikes_bar(a,tp,2))*1.3, 'r*')
       end
   end
end
%Firing rates during natural ROC
plot([0.5 1.5], [3.4555 3.4555], 'b--'); %PFC
plot([1.5 2.5], [3.1227 3.1227], 'r--'); %8A
plot([2.5 3.5], [3.5027 3.5027], 'g--'); %PPC
plot([3.5 4.5], [5.0097 5.0097], 'm--'); %STG
ylim([0.5 7])
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'})
ylabel('Firing Rate (spikes/sec)')
title('High Current - Nonparametric')
subplot(2,1,2);
bar(squeeze(spikes_bar(:,:,1))); hold on;
%pre-stim
NEG = squeeze(spikes_bar(:,1,1)-spikes_CI_L(:,1,1)); POS = squeeze(spikes_CI_H(:,1,1)-squeeze(spikes_bar(:,1,1)));
errorbar([1:4]-0.31, squeeze(spikes_bar(:,1,1)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)
%stim
NEG = squeeze(spikes_bar(:,2,1)-spikes_CI_L(:,2,1)); POS = squeeze(spikes_CI_H(:,1,1)-squeeze(spikes_bar(:,1,1)));
errorbar([1:4]-0.16, squeeze(spikes_bar(:,2,1)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)
%post1
NEG = squeeze(spikes_bar(:,3,1)-spikes_CI_L(:,3,1)); POS = squeeze(spikes_CI_H(:,1,1)-squeeze(spikes_bar(:,1,1)));
errorbar([1:4], squeeze(spikes_bar(:,3,1)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)
%post2
NEG = squeeze(spikes_bar(:,4,1)-spikes_CI_L(:,4,1)); POS = squeeze(spikes_CI_H(:,1,1)-squeeze(spikes_bar(:,1,1)));
errorbar([1:4]+0.16, squeeze(spikes_bar(:,4,1)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)
%post3
NEG = squeeze(spikes_bar(:,5,1)-spikes_CI_L(:,5,1)); POS = squeeze(spikes_CI_H(:,1,1)-squeeze(spikes_bar(:,1,1)));
errorbar([1:4]+0.31, squeeze(spikes_bar(:,5,1)), NEG, POS, 'k', 'linestyle', 'none', 'linewidth', 1.5)

for a = 1:4
   for tp = 2:5
       if spikes_stats(a,tp,1)<0.01
           plot(a+xpos(tp), squeeze(spikes_bar(a,tp,1))*1.3, 'r*')
       end
   end
end
%Firing rates during natural ROC
plot([0.5 1.5], [3.4555 3.4555], 'b--'); %PFC
plot([1.5 2.5], [3.1227 3.1227], 'r--'); %8A
plot([2.5 3.5], [3.5027 3.5027], 'g--'); %PPC
plot([3.5 4.5], [5.0097 5.0097], 'm--'); %STG
ylim([0.5 7])
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'})
ylabel('Firing Rate (spikes/sec)')
title('Low Current - Nonparametric')


t1 = nearest(time_axis, 1);
t2 = nearest(time_axis, 28);

[h p] = signtest(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2)-5.0097)
[h p] = signtest(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2)-3.5)
[h p] = signtest(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2)-3.12)
[h p] = signtest(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2)-3.45)


[h p] = ttest(nanmean(spike_rates_bs_sessions{4,2}(:,t1:t2),2)-5.0097)


[h p] = ttest(nanmean(spike_rates_bs_sessions{3,2}(:,t1:t2),2)-3.5)
[h p] = ttest(nanmean(spike_rates_bs_sessions{2,2}(:,t1:t2),2)-3.12)
[h p] = ttest(nanmean(spike_rates_bs_sessions{1,2}(:,t1:t2),2)-3.45)
