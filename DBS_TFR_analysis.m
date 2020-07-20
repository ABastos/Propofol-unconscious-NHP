
cd \\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolWakeUp

filestoload = ls('\\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolWakeUp\*.mat');


filestoload(end,:) = [];

% for s = 1:11 %size(filestoload,1)
for s = 1:size(filestoload,1)
    cd('\\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolWakeUp\')
    
    load(deblank(filestoload(s,:)), 'trialInfo', 'sessionInfo', 'electrodeInfo', 'lfp', 'lfpSchema')
    
    
    isStandardAllAmps = trialInfo.dbs_isWakeUpTest & ...      % Is it a wake-up test?
        ismember(trialInfo.dbs_testType,{'WUT7','WUT7ab','WUT7abU'}) & ... % WUT7/ab family of tests
        ... (trialInfo.dbs_amplitude > 0.5) & ... % Amplitude > 0.5 mAmp
        (trialInfo.dbs_frequency == 180) & ...% Frequency = 180 hz
        (trialInfo.dbs_duration > 25) & ...   % 20 s < Duration < 35 s
        (trialInfo.dbs_duration < 35) & ...
        (trialInfo.dbs_previousITI > 120) ; %previous stim was greater than 2 minutes ago
    
    standardTrialIndx = find(isStandardAllAmps);
    
    numtr = length(standardTrialIndx);
    ntime = size(lfp,1);
    nchannels = size(lfpSchema.index{2},2); %size(electrodeInfo.electrode,1); %+size(ain,2);
    
    %find correspondence between lfp field and electrodeInfo:
    chnl_mapping = zeros(nchannels,1);
    for c = 1:nchannels
        chnl_mapping(c) = find(ismember(electrodeInfo.chnlID, lfpSchema.index{2}{c}));
    end
    
    %get data into DBS trial format
    data =[];
    data.trialinfo = zeros(numtr,4);
    for tr = 1:numtr
        data.time{tr} = -30:0.001:125;
    
        data.trial{tr} = zeros(nchannels, length(data.time{1})); %num channel x num time
        
        tr_start = trialInfo.dbs_stimOn(standardTrialIndx(tr))*1000;
        
        data.trial{tr}(:,:) = squeeze(lfp(tr_start-(30*1000):tr_start+(155*1000),:))';
        
        data.trialinfo(tr,1) = standardTrialIndx(tr);
        data.trialinfo(tr,2) = trialInfo.dbs_amplitude(standardTrialIndx(tr));
        data.trialinfo(tr,3) = trialInfo.dbs_wakeScore(standardTrialIndx(tr));
        data.trialinfo(tr,4) = trialInfo.dbs_duration(standardTrialIndx(tr));
    end
    
    if strcmp(sessionInfo.subject, 'Mary')
        data.trialinfo(:,5) = 2.*ones(numtr,1);
    else
        data.trialinfo(:,5) = ones(numtr,1);
    end
    
    data.fsample = 1000;
    data.label = {};
    %make each channel unique
    for c = 1:nchannels
        data.label{c} = [electrodeInfo.area{chnl_mapping(c)} int2str(electrodeInfo.channel(chnl_mapping(c)))];
    end
    

%     cfg = [];
%     cfg.viewmode = 'vertical';
% %     cfg.channel = [1 size(electrodeInfo.electrode,1)+1:length(data.label)];
%     cfg.channel = [1 65 126 190];
%     % cfg.chanscale = stdchan([1 size(electrodeInfo.electrode,1)+1:length(data.label)]);
%     
%     ft_databrowser(cfg, data)    
    
    

    electrodeInfo.electrode = electrodeInfo.electrode(chnl_mapping);
    electrodeInfo.chnlID = electrodeInfo.chnlID(chnl_mapping);
    electrodeInfo.channel = electrodeInfo.channel(chnl_mapping);
    electrodeInfo.session = electrodeInfo.session(chnl_mapping);
    electrodeInfo.file = electrodeInfo.file(chnl_mapping);
    electrodeInfo.hasLFP = electrodeInfo.hasLFP(chnl_mapping);
    electrodeInfo.numUnits = electrodeInfo.numUnits(chnl_mapping);
    electrodeInfo.hemisphere = electrodeInfo.hemisphere(chnl_mapping);
    electrodeInfo.area = electrodeInfo.area(chnl_mapping);
    electrodeInfo.array = electrodeInfo.array(chnl_mapping);
    electrodeInfo.NSP = electrodeInfo.NSP(chnl_mapping);
    electrodeInfo.gridLoc = electrodeInfo.gridLoc(chnl_mapping,:);
    
%     [data_bipolar2 ] = derive_400um_bipolar_Utah(data, electrodeInfo)
    

%     cfg = [];
%     cfg.viewmode = 'vertical';
%     cfg.channel = [1 45 75 100];
%     
%     ft_databrowser(cfg, data_bipolar2)    
    
    
    %remove the ERP from each trial condition (by amount of current)
    
    current_levels = [0.5 1.0 1.5 2.0 2.5];

    for curr = 1:length(current_levels)
        current_trials = find(data.trialinfo(:,2) == current_levels(curr));
        if ~isempty(current_trials)
            cfg = [];
            cfg.trials = current_trials;

            tmp = ft_selectdata(cfg, data);
            tmp = ft_timelockanalysis([], tmp);

            for tr = 1:length(current_trials)
                data.trial{current_trials(tr)} = data.trial{current_trials(tr)} - tmp.avg;
            end
        end
    end    
    
%     %If doing bipolar:
%     data_bipolar3 = data_bipolar2;
%     for curr = 1:length(current_levels)
%         
%         current_trials = find(data_bipolar2.trialinfo(:,2) == current_levels(curr));
%         if ~isempty(current_trials)
%             cfg = [];
%             cfg.trials = current_trials;
% 
%             tmp = ft_selectdata(cfg, data_bipolar2);
% 
%             tmp = ft_timelockanalysis([], tmp);
% 
%             for tr = 1:length(current_trials)
%                 data_bipolar3.trial{current_trials(tr)} = data_bipolar3.trial{current_trials(tr)} - tmp.avg;
%             end
%         end
%     end
    
%     cfg = [];
%     cfg.viewmode = 'vertical';
%     cfg.channel = [1 45 75 100];
%     
%     ft_databrowser(cfg, data_bipolar3) 
    
    PFC_indx = [];
    FEF_indx = [];
    PPC_indx = [];
    STG_indx = [];
    for c = 1:length(data.label)
        if strcmp(data.label{c}(1:2), '7b')
            PPC_indx = [PPC_indx; c];
        elseif strcmp(data.label{c}(1:3), 'FEF')
            FEF_indx = [FEF_indx; c];
        elseif strcmp(data.label{c}(1:3), 'CPB')
            STG_indx = [STG_indx; c];
        elseif strcmp(data.label{c}(1:5), 'vlPFC')
            PFC_indx = [PFC_indx; c];
        end
    end
    
    
    foi = logspace(-0.75, 2.3, 50);
    toi = -25:0.5:150;
    ntimefreq = length(toi);
    
    ALL_TFR = {};
    ALL_TFRbs = {};
    
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.toi = toi;
    cfg.t_ftimwin = 5.*ones(length(foi),1);
    cfg.channel = PFC_indx;
    cfg.taper = 'hanning';
    % cfg.tapsmofrq = 6;
    cfg.output = 'pow';
    cfg.keeptrials = 'yes';
    cfg.foi = foi;
    ALL_TFR{1} = ft_freqanalysis(cfg, data)
    cfg.channel = FEF_indx;
    ALL_TFR{2} = ft_freqanalysis(cfg, data)
    cfg.channel = PPC_indx;
    ALL_TFR{3} = ft_freqanalysis(cfg, data)
    cfg.channel = STG_indx;
    ALL_TFR{4} = ft_freqanalysis(cfg, data)    
%     
    
    cfg = [];
    cfg.baseline = [-25 -3];
    cfg.baselinetype = 'db';
    ALL_TFRbs{1} = ft_freqbaseline(cfg, ALL_TFR{1})
    ALL_TFRbs{2} = ft_freqbaseline(cfg, ALL_TFR{2})
    ALL_TFRbs{3} = ft_freqbaseline(cfg, ALL_TFR{3})
    ALL_TFRbs{4} = ft_freqbaseline(cfg, ALL_TFR{4})
%     
    

    %create time axis
    time_axis_baseline = ALL_TFRbs{1}.time;
    xmarkers = 1:10:length(time_axis_baseline);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis_baseline(xmarkers(indx)));
        indx=indx+1;
    end
    
    
    %create frequency axis
    freqaxis = ALL_TFRbs{1}.freq;
    ymarkers = 1:5:length(freqaxis);
    ylabels = {};
    indx=1;
    for y = ymarkers
        ylabels{indx} = num2str(freqaxis(ymarkers(indx)), '%0.3g');
        indx=indx+1;
    end
%     
%     for c = 1:length(ALL_TFRbs{1}.label)
% 
%         %Hi vs low current
%         current_levels = [0.5 1.0 1.5 2.0 2.5];
%         current_trials = find(ALL_TFRbs{1}.trialinfo(:,2) == current_levels(1));
%         figure; 
%         subplot(2,1,1);
%         imagesc(squeeze(mean(ALL_TFRbs{1}.powspctrm(current_trials,c,:,:),1))); hold on;
%         set(gca, 'ydir', 'normal')
%         colormap('Jet')
%         set(gca, 'xtick', xmarkers)
%         set(gca, 'xticklabel', xlabels)
%         xlabel('Time since Drug start (minutes)')
%         set(gca, 'ytick', ymarkers)
%         set(gca, 'yticklabel', ylabels)
%         ylabel('Frequency (Hz)')
%         xline(nearest(ALL_TFRbs{1}.time, 0), 'k', 'LineWidth', 2)
%         dur = mean(ALL_TFRbs{1}.trialinfo(:,4));
%         durbin = nearest(ALL_TFRbs{1}.time, dur);
%         xline(durbin, 'k', 'LineWidth', 2)
%         caxis([-5 5])
%         colorbar
%         title(['PFC - low current, chan = ', int2str(c)])
% 
%         subplot(2,1,2);
%         current_trials = find(ALL_TFRbs{1}.trialinfo(:,2) == current_levels(5));
%         imagesc(squeeze(mean(ALL_TFRbs{1}.powspctrm(current_trials,c,:,:),1))); hold on;
%         set(gca, 'ydir', 'normal')
%         colormap('Jet')
%         set(gca, 'xtick', xmarkers)
%         set(gca, 'xticklabel', xlabels)
%         xlabel('Time since Drug start (minutes)')
%         set(gca, 'ytick', ymarkers)
%         set(gca, 'yticklabel', ylabels)
%         ylabel('Frequency (Hz)')
%         xline(nearest(ALL_TFRbs{1}.time, 0), 'k', 'LineWidth', 2)
%         dur = mean(ALL_TFRbs{1}.trialinfo(:,4));
%         durbin = nearest(ALL_TFRbs{1}.time, dur);
%         xline(durbin, 'k', 'LineWidth', 2)
%         caxis([-5 5])
%         colorbar
%         title(['PFC - high currentm chan = ',  int2str(c)] )
% 
%         
%         pause
%         close all
%         
%     end
    

    if s >= 1 & s <= 11
        %mary
        low_current = [0.5 1.0];
        high_current = [1.5 2.5];
%         tmp_trialinfo(:,5) = 2.*ones(size(ALL_TFRbs{1}.trialinfo(:,:),1),1);
    else
        %jones
        low_current = [0.5 NaN];
        high_current = [1.0 1.5];
%         tmp_trialinfo(:,5) = ones(size(ALL_TFRbs{1}.trialinfo(:,:),1),1);
    end
    
    %Hi vs low current
%     current_levels = [0.5 1.0 1.5 2.0 2.5];
    current_trials = find(ALL_TFRbs{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs{1}.trialinfo(:,2) == low_current(2));
    figure; 
    subplot(2,1,1);
    imagesc(squeeze(mean(mean(ALL_TFRbs{1}.powspctrm(current_trials,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-4 4])
    colorbar
    title(['PFC - low current, All chans'])

    current_trials = find(ALL_TFRbs{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs{1}.trialinfo(:,2) == high_current(2));
    subplot(2,1,2);
    imagesc(squeeze(mean(mean(ALL_TFRbs{1}.powspctrm(current_trials,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-4 4])
    colorbar
    title(['PFC - high current, All chans'] )


    cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')
    save([sessionInfo.session '_DBS_TFR_2minspost'], 'ALL_TFR', 'ALL_TFRbs', 'sessionInfo','trialInfo', 'standardTrialIndx',  '-v7.3')
    
    clear data lfp ALL_TFR ALL_TFRbs sessionInfo data_bipolar2 data_bipolar3
end

return
%%
% filestoload = ls('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\*DBS_TFR*.mat');
filestoload = ls('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\*DBS_TFR_2minspost_uni*.mat')

current_levels = [0.5 1.0 1.5 2.0 2.5];

ALL_TFRbs_sessions = {};
ALL_TFR_sessions = {};
for s = 1:size(filestoload,1)

    load(deblank(filestoload(s,:)), 'ALL_TFR')
    
%     tmp_trialinfo = zeros(size(ALL_TFRbs{1}.trialinfo(:,:),1), 5);
%     tmp_trialinfo(:,1:4) = ALL_TFRbs{1}.trialinfo;
    

    if s >= 1 & s <= 11
        %mary
        low_current = [0.5 1.0];
        high_current = [1.5 2.5];
%         tmp_trialinfo(:,5) = 2.*ones(size(ALL_TFRbs{1}.trialinfo(:,:),1),1);
    else
        %jones
        low_current = [0.5 NaN];
        high_current = [1.0 1.5];
%         tmp_trialinfo(:,5) = ones(size(ALL_TFRbs{1}.trialinfo(:,:),1),1);
    end
%     ALL_TFRbs{1}.trialinfo = tmp_trialinfo;
%     ALL_TFRbs{2}.trialinfo = tmp_trialinfo;
%     ALL_TFRbs{3}.trialinfo = tmp_trialinfo;
%     ALL_TFRbs{4}.trialinfo = tmp_trialinfo;
    


    cfg = [];
    cfg.baseline = [-10 -3];
    cfg.baselinetype = 'db';
    ALL_TFRbs{1} = ft_freqbaseline(cfg, ALL_TFR{1})
    ALL_TFRbs{2} = ft_freqbaseline(cfg, ALL_TFR{2})
    ALL_TFRbs{3} = ft_freqbaseline(cfg, ALL_TFR{3})
    ALL_TFRbs{4} = ft_freqbaseline(cfg, ALL_TFR{4})
    
    areas = {'PFC', '8A', 'PPC', 'STG'};
    if s ==1
        for a = 1:4
            ALL_TFRbs{a}.powspctrm = mean(ALL_TFRbs{a}.powspctrm,2);
            ALL_TFRbs{a}.label = areas(a);
            ALL_TFRbs_sessions{a} = ALL_TFRbs{a};
            
            ALL_TFR{a}.powspctrm = mean(ALL_TFR{a}.powspctrm,2);
            ALL_TFR{a}.label = areas(a);            
            ALL_TFR_sessions{a} = ALL_TFR{a};
        end
        
    else
        for a = 1:4
            ALL_TFRbs{a}.powspctrm = mean(ALL_TFRbs{a}.powspctrm,2);
            ALL_TFRbs{a}.label = areas(a);            
            ALL_TFRbs_sessions{a} = ft_appendfreq([], ALL_TFRbs_sessions{a}, ALL_TFRbs{a});
            
            ALL_TFR{a}.powspctrm = mean(ALL_TFR{a}.powspctrm,2);
            ALL_TFR{a}.label = areas(a);            
            ALL_TFR_sessions{a} = ft_appendfreq([], ALL_TFR_sessions{a}, ALL_TFR{a});            
        end
    end
    
    %create time axis
    time_axis_baseline = ALL_TFRbs{1}.time;
    xmarkers = 1:10:length(time_axis_baseline);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis_baseline(xmarkers(indx)));
        indx=indx+1;
    end
    
    
    %create frequency axis
    freqaxis = ALL_TFRbs{1}.freq;
    ymarkers = 1:5:length(freqaxis);
    ylabels = {};
    indx=1;
    for y = ymarkers
        ylabels{indx} = num2str(freqaxis(ymarkers(indx)), '%0.3g');
        indx=indx+1;
    end    
    
    %Hi vs low current
    xlim1 = nearest(time_axis_baseline, -10);
    xlim2 = nearest(time_axis_baseline, 135);
    
    current_trials = find(ALL_TFRbs{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs{1}.trialinfo(:,2) == low_current(2));
    figure; 
    subplot(2,1,1);
    imagesc(squeeze(mean(mean(ALL_TFRbs{2}.powspctrm(current_trials,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-5 5])
    xlim([xlim1 xlim2])
    colorbar
    title(['PFC - low current, All chans, N = ' int2str(length(current_trials)) ' session: ' int2str(s)])

    subplot(2,1,2);
    current_trials = find(ALL_TFRbs{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs{1}.trialinfo(:,2) == high_current(2));
    imagesc(squeeze(mean(mean(ALL_TFRbs{2}.powspctrm(current_trials,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-5 5])
    xlim([xlim1 xlim2])
    colorbar
    title(['PFC - high current, All chans, N = ' int2str(length(current_trials)) ' session: ' int2str(s)])

end


%%

    %create time axis
    time_axis_baseline = ALL_TFRbs{1}.time;
    xlim1 = nearest(time_axis_baseline, -10);
    xlim2 = nearest(time_axis_baseline, 135);
    
    if s >= 1 & s <= 11
        %mary
        low_current = [0.5 1.0];
        high_current = [1.5 2.5];
%         tmp_trialinfo(:,5) = 2.*ones(size(ALL_TFRbs{1}.trialinfo(:,:),1),1);
    else
        %jones
        low_current = [0.5 NaN];
        high_current = [1.0 1.5];
%         tmp_trialinfo(:,5) = ones(size(ALL_TFRbs{1}.trialinfo(:,:),1),1);
    end
    
    
        low_current = [0.5 1.0];
        high_current = [1.5 2.5];
%MAry
    current_trials = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 2);
    figure; 
    subplot(2,1,1);
    imagesc(squeeze(mean(mean(ALL_TFRbs_sessions{1}.powspctrm(current_trials,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs_sessions{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs_sessions{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs_sessions{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-5 5])
    xlim([xlim1 xlim2])
    colorbar
    title(['PFC (Mary) - low current, All chans, N = ' int2str(length(current_trials)) ])

    subplot(2,1,2);
    current_trials = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 2);
    imagesc(squeeze(mean(nanmean(ALL_TFRbs_sessions{1}.powspctrm(current_trials,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs_sessions{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs_sessions{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs_sessions{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-5 5])
    xlim([xlim1 xlim2])
    colorbar
    title(['PFC (Mary) - high current, All chans, N = ' int2str(length(current_trials)) ])

  
%Jones
        low_current = [0.5 NaN];
        high_current = [1.0 1.5];
        
    current_trials = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 1);
    figure; 
    subplot(2,1,1);
    imagesc(squeeze(nanmean(nanmean(ALL_TFRbs_sessions{1}.powspctrm(current_trials,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs_sessions{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs_sessions{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs_sessions{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-5 5])
    xlim([xlim1 xlim2])
    colorbar
    title(['PFC (Jones) - low current, All chans, N = ' int2str(length(current_trials)) ])

    subplot(2,1,2);
    current_trials = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 1);
    imagesc(squeeze(nanmean(nanmean(ALL_TFRbs_sessions{1}.powspctrm(current_trials,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs_sessions{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs_sessions{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs_sessions{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-5 5])
    xlim([xlim1 xlim2])
    colorbar
    title(['PFC (Jones) - high current, All chans, N = ' int2str(length(current_trials)) ])



%Combined monkeys
        low_current = [0.5 1.0];
        %MAry
    current_trials1 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 2);

%Jones
        low_current = [0.5 NaN];        
    current_trials2 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 1);
    
current_trials_low = [current_trials1; current_trials2];
%%
%Combined monkeys
        high_current = [1.5 2.5];
%MAry
    current_trials1 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 2);

%Jones
        high_current = [1.0 1.5];
        
    current_trials2 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 1);
    
current_trials_high = [current_trials1; current_trials2];
    a=4
    area_list = {'PFC', '8A', 'PPC', 'STG'};
    figure; 
    subplot(2,1,1);
    imagesc(squeeze(nanmean(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_low,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs_sessions{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs_sessions{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs_sessions{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-5 5])
    xlim([xlim1 xlim2])
    colorbar
    title([area_list{a} ' (Combined) - low current, All chans, N = ' int2str(length(current_trials_low)) ])

    subplot(2,1,2);
    imagesc(squeeze(nanmean(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_high,:,:,:),1),2))); hold on;
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    ylabel('Frequency (Hz)')
    xline(nearest(ALL_TFRbs_sessions{1}.time, 0), 'k', 'LineWidth', 2)
    dur = mean(ALL_TFRbs_sessions{1}.trialinfo(:,4));
    durbin = nearest(ALL_TFRbs_sessions{1}.time, dur);
    xline(durbin, 'k', 'LineWidth', 2)
    caxis([-5 5])
    xlim([xlim1 xlim2])
    colorbar
    title([area_list{a}  ' (Combined) - high current, All chans, N = ' int2str(length(current_trials_high)) ])
    
%%
%Bar plots
TFR_bar = zeros(4,4,3,2); %Dim is area x time period ( stim, post1/2/3) x freqband x current (low,high)

TFR_bar_sem = zeros(4,4,2);
bs1 = nearest(ALL_TFRbs_sessions{1}.time, -30);
bs2 = nearest(ALL_TFRbs_sessions{1}.time, -3);
for freqband = 1:3
    if freqband == 1
        %SF
        f1 = nearest(ALL_TFRbs{1}.freq, 0.1);
        f2 = nearest(ALL_TFRbs{1}.freq, 1.5);
    elseif freqband == 2
        %Beta
        f1 = nearest(ALL_TFRbs{1}.freq, 10);
        f2 = nearest(ALL_TFRbs{1}.freq, 25);
    elseif freqband == 3
        %Gamma
        f1 = nearest(ALL_TFRbs{1}.freq, 50);
        f2 = nearest(ALL_TFRbs{1}.freq, 200);
    end
    for timeperiod = 1:4
        if timeperiod == 1
            tbin1 = nearest(ALL_TFRbs_sessions{1}.time, 3);
            tbin2 = nearest(ALL_TFRbs_sessions{1}.time, 25);
        elseif timeperiod == 2
            tbin1 = nearest(ALL_TFRbs_sessions{1}.time, 32);
            tbin2 = nearest(ALL_TFRbs_sessions{1}.time, 60);
        elseif timeperiod == 3
            tbin1 = nearest(ALL_TFRbs_sessions{1}.time, 61);
            tbin2 = nearest(ALL_TFRbs_sessions{1}.time, 90);
        elseif timeperiod == 4
            tbin1 = nearest(ALL_TFRbs_sessions{1}.time, 91);
            tbin2 = nearest(ALL_TFRbs_sessions{1}.time, 120);
        end
        
        for a = 1:4
            %Low Current
            ntrials = length(current_trials_low);
            tmp_pow1 = squeeze(nanmean(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_low,:,f1:f2,tbin1:tbin2),3),4));
            bs = squeeze(nanmean(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_low,:,f1:f2,bs1:bs2),3),4));
            tmp_pow1 = tmp_pow1 - bs;
            TFR_bar(a,timeperiod,freqband,1) = nanmean(tmp_pow1);
            TFR_bar_sem(a,timeperiod,freqband,1) = nanstd(tmp_pow1)./sqrt(ntrials);

            %High Current
            ntrials = length(current_trials_high);   
            tmp_pow2 = squeeze(nanmean(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_high,:,f1:f2,tbin1:tbin2),3),4));
            bs = squeeze(nanmean(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_high,:,f1:f2,bs1:bs2),3),4));
            tmp_pow2 = tmp_pow2 - bs;            
            TFR_bar(a,timeperiod,freqband,2) = nanmean(tmp_pow2);
            TFR_bar_sem(a,timeperiod,freqband,2) = nanstd(tmp_pow2)./sqrt(ntrials);
            
        end
        
    end
end


% TFR_bar = zeros(4,4,3,2); %Dim is area x time period ( stim, post1/2/3) x freqband x current (low,high)
%8A, post3, SF, low
TFR_bar(2,3,1,1)


%PFC, post3, SF, high
TFR_bar(1,4,1,2)
%8A, post3, SF, high
TFR_bar(2,4,1,2)
%PPC, post3, SF, high
TFR_bar(3,4,1,2)
%STG, post3, SF, high
TFR_bar(4,4,1,2)

%PFC, post3, SF, low
TFR_bar(1,4,1,1)
%8A, post3, SF, low
TFR_bar(2,4,1,1)
%PPC, post3, SF, low
TFR_bar(3,4,1,1)
%STG, post3, SF, low
TFR_bar(4,4,1,1)



%PFC, dur, Beta/Gamma, high
TFR_bar(1,1,3,2)
%8A, dur, Beta/Gamma, high
TFR_bar(2,1,3,2)
%PPC, dur, Beta/Gamma, high
TFR_bar(3,1,3,2)
%STG, dur, Beta/Gamma, high
TFR_bar(4,1,3,2)

%PFC, dur, Beta/Gamma, low
TFR_bar(1,4,2,1)
%8A, dur, Beta/Gamma, low
TFR_bar(2,4,2,1)
%PPC, dur, Beta/Gamma, low
TFR_bar(3,4,2,1)
%STG, dur, Beta/Gamma, low
TFR_bar(4,4,2,1)



figure;
subplot(2,1,1);
bar([1:4], squeeze(TFR_bar(:,:,1,2))); hold on
errorbar([1:4]-0.27, squeeze(TFR_bar(:,1,1,2)), 2.*squeeze(TFR_bar_sem(:,1,1,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]-0.09, squeeze(TFR_bar(:,2,1,2)), 2.*squeeze(TFR_bar_sem(:,2,1,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.09, squeeze(TFR_bar(:,3,1,2)), 2.*squeeze(TFR_bar_sem(:,3,1,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.27, squeeze(TFR_bar(:,4,1,2)), 2.*squeeze(TFR_bar_sem(:,4,1,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
% ylim([-8 0.5])
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'})
ylabel('dB Change from baseline')
title('High Current - SF')
subplot(2,1,2);
bar([1:4], squeeze(TFR_bar(:,:,1,1))); hold on
errorbar([1:4]-0.27, squeeze(TFR_bar(:,1,1,1)), 2.*squeeze(TFR_bar_sem(:,1,1,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]-0.09, squeeze(TFR_bar(:,2,1,1)), 2.*squeeze(TFR_bar_sem(:,2,1,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.09, squeeze(TFR_bar(:,3,1,1)), 2.*squeeze(TFR_bar_sem(:,3,1,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.27, squeeze(TFR_bar(:,4,1,1)), 2.*squeeze(TFR_bar_sem(:,4,1,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
% ylim([-8 0.5])
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'})
ylabel('dB Change from baseline')
title('Low Current - SF')



figure;
subplot(2,1,1);
bar([1:4], squeeze(TFR_bar(:,:,3,2))); hold on
errorbar([1:4]-0.27, squeeze(TFR_bar(:,1,3,2)), 2.*squeeze(TFR_bar_sem(:,1,3,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]-0.09, squeeze(TFR_bar(:,2,3,2)), 2.*squeeze(TFR_bar_sem(:,2,3,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.09, squeeze(TFR_bar(:,3,3,2)), 2.*squeeze(TFR_bar_sem(:,3,3,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.27, squeeze(TFR_bar(:,4,3,2)), 2.*squeeze(TFR_bar_sem(:,4,3,2)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
% ylim([-0.5 6])
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'})
ylabel('dB Change from baseline')
title('High Current - Gamma')
subplot(2,1,2);
bar([1:4], squeeze(TFR_bar(:,:,3,1))); hold on
errorbar([1:4]-0.27, squeeze(TFR_bar(:,1,3,1)), 2.*squeeze(TFR_bar_sem(:,1,3,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]-0.09, squeeze(TFR_bar(:,2,3,1)), 2.*squeeze(TFR_bar_sem(:,2,3,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.09, squeeze(TFR_bar(:,3,3,1)), 2.*squeeze(TFR_bar_sem(:,3,3,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
errorbar([1:4]+0.27, squeeze(TFR_bar(:,4,3,1)), 2.*squeeze(TFR_bar_sem(:,4,3,1)), 'k', 'linestyle', 'none', 'linewidth', 1.5)
% ylim([-2.5 5.1])
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'})
ylabel('dB Change from baseline')
title('Low Current - Gamma')


%%
%Change from baseline spectra


%Combined monkeys
        low_current = [0.5 1.0];
        %MAry
    current_trials1 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 2);

%Jones
        low_current = [0.5 NaN];        
    current_trials2 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 1);
    
current_trials_low = [current_trials1; current_trials2];

%Combined monkeys
        high_current = [1.5 2.5];
%MAry
    current_trials1 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 2);

%Jones
        high_current = [1.0 1.5];
        
    current_trials2 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 1);
    
current_trials_high = [current_trials1; current_trials2];

%defined behaviorally:
% current_trials_high = (find(ALL_TFR_sessions{1}.trialinfo(:,3)>=4 & ALL_TFR_sessions{1}.trialinfo(:,3)<=6))
% current_trials_low = (find(ALL_TFR_sessions{1}.trialinfo(:,3)>=0 & ALL_TFR_sessions{1}.trialinfo(:,3)<=2))


TFR_spectra = zeros(4,4,50,2); %Dim is area x time period ( stim, post1/2/3) x freqband x current (low,high)

TFR_spectra_sem = zeros(4,4,50,2);
bs1 = nearest(ALL_TFRbs_sessions{1}.time, -30);
bs2 = nearest(ALL_TFRbs_sessions{1}.time, -3);

for timeperiod = 1:4
    if timeperiod == 1
        tbin1 = nearest(ALL_TFRbs_sessions{1}.time, 3);
        tbin2 = nearest(ALL_TFRbs_sessions{1}.time, 25);%         tmp_pow2 = squeeze(nanmean(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_high,:,f1:f2,tbin1:tbin2),3),4));
%         bs = squeeze(nanmean(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_high,:,f1:f2,bs1:bs2),3),4));
%         tmp_pow2 = tmp_pow2 - bs;

    elseif timeperiod == 2
        tbin1 = nearest(ALL_TFRbs_sessions{1}.time, 32);
        tbin2 = nearest(ALL_TFRbs_sessions{1}.time, 60);
    elseif timeperiod == 3
        tbin1 = nearest(ALL_TFRbs_sessions{1}.time, 61);
        tbin2 = nearest(ALL_TFRbs_sessions{1}.time, 90);
    elseif timeperiod == 4
        tbin1 = nearest(ALL_TFRbs_sessions{1}.time, 91);
        tbin2 = nearest(ALL_TFRbs_sessions{1}.time, 120);
    end
    
    for a = 1:4
        %Low Current
        ntrials = length(current_trials_low);
        tmp_pow1 = squeeze(nanmean(ALL_TFR_sessions{a}.powspctrm(current_trials_low,:,:,tbin1:tbin2),4));
        bs = squeeze(nanmean(ALL_TFR_sessions{a}.powspctrm(current_trials_low,:,:,bs1:bs2),4));
        tmp_pow1 = 10.*log10(tmp_pow1./bs);
        TFR_spectra(a,timeperiod,:,1) = nanmean(tmp_pow1);
        TFR_spectra_sem(a,timeperiod,:,1) = nanstd(tmp_pow1)./sqrt(ntrials);
        
        %High Current
%         tmp_pow2 = squeeze(nanmean(ALL_TFRbs_sessions{a}.powspctrm(current_trials_high,:,:,tbin1:tbin2),4));
        tmp_pow2 = squeeze(nanmean(ALL_TFR_sessions{a}.powspctrm(current_trials_high,:,:,tbin1:tbin2),4));
        bs = squeeze(nanmean(ALL_TFR_sessions{a}.powspctrm(current_trials_high,:,:,bs1:bs2),4));
        tmp_pow2 = 10.*log10(tmp_pow2./bs);
        ntrials = length(current_trials_high);
        TFR_spectra(a,timeperiod,:,2) = nanmean(tmp_pow2);
        TFR_spectra_sem(a,timeperiod,:,2) = nanstd(tmp_pow2)./sqrt(ntrials);
        
    end
    
end


figure; shadedErrorBar(freqaxis, squeeze(TFR_spectra(1,1,:,1)), 2.*squeeze(TFR_spectra_sem(1,1,:,1)), 'b');
hold on;
shadedErrorBar(freqaxis, squeeze(TFR_spectra(1,1,:,2)), 2.*squeeze(TFR_spectra_sem(1,1,:,2)), 'r');
set(gca, 'xscale', 'log')
yline(0, 'k--')
xlim([freqaxis(1) freqaxis(end)])
title('PFC')

timelabel = {'stim','post1','post2','post3'};
colortoplot = {'r', 'y', 'k', 'b'};
figure; 
for t = 1:4
% figure; 
shadedErrorBar(freqaxis, squeeze(TFR_spectra(4,t,:,2)), 2.*squeeze(TFR_spectra_sem(4,t,:,2)), colortoplot{t});
hold on;
% shadedErrorBar(freqaxis, squeeze(TFR_spectra(2,t,:,2)), 2.*squeeze(TFR_spectra_sem(2,t,:,2)), 'r');
set(gca, 'xscale', 'log')
yline(0, 'k--')
xlim([freqaxis(1) freqaxis(end)])
title(['PFC, time period: ' timelabel{t}])
end



timelabel = {'stim','post1','post2','post3'};
colortoplot = {'b', 'r', 'g', 'm'};
figure; 
for a = [3 4 1 2]
% figure; 
shadedErrorBar(freqaxis, squeeze(TFR_spectra(a,1,:,2)), 2.*squeeze(TFR_spectra_sem(a,1,:,2)), colortoplot{a});
hold on;
% shadedErrorBar(freqaxis, squeeze(TFR_spectra(2,t,:,2)), 2.*squeeze(TFR_spectra_sem(2,t,:,2)), 'r');
set(gca, 'xscale', 'log')
yline(0, 'k--')
xlim([freqaxis(1) freqaxis(end)])
title(['stim period (High Current)'])
ylim([-10 10])
xlabel('Frequency (Hz)')
ylabel('dB Change from pre-stim baseline')
end


timelabel = {'stim','post1','post2','post3'};
colortoplot = {'b', 'r', 'g', 'm'};
figure; 
for a = [3 4 1 2]
% figure; 
shadedErrorBar(freqaxis, squeeze(TFR_spectra(a,2,:,2)), 2.*squeeze(TFR_spectra_sem(a,2,:,2)), colortoplot{a});
hold on;
% shadedErrorBar(freqaxis, squeeze(TFR_spectra(2,t,:,2)), 2.*squeeze(TFR_spectra_sem(2,t,:,2)), 'r');
set(gca, 'xscale', 'log')
yline(0, 'k--')
xlim([freqaxis(1) freqaxis(end)])
title([ 'post1 period (High Current)'])
ylim([-10 10])
xlabel('Frequency (Hz)')
ylabel('dB Change from pre-stim baseline')
end

figure; shadedErrorBar(freqaxis, squeeze(TFR_spectra(3,1,:,1)), 2.*squeeze(TFR_spectra_sem(3,1,:,1)), 'b');
hold on;
shadedErrorBar(freqaxis, squeeze(TFR_spectra(3,1,:,2)), 2.*squeeze(TFR_spectra_sem(3,1,:,2)), 'r');
set(gca, 'xscale', 'log')
yline(0, 'k--')
xlim([freqaxis(1) freqaxis(end)])
title('PPC')

figure; shadedErrorBar(freqaxis, squeeze(TFR_spectra(4,1,:,1)), 2.*squeeze(TFR_spectra_sem(4,1,:,1)), 'b');
hold on;
shadedErrorBar(freqaxis, squeeze(TFR_spectra(4,1,:,2)), 2.*squeeze(TFR_spectra_sem(4,1,:,2)), 'r');
set(gca, 'xscale', 'log')
yline(0, 'k--')
xlim([freqaxis(1) freqaxis(end)])
title('STG')
%All areas:
figure; 
shadedErrorBar(freqaxis, squeeze(TFR_spectra(1,1,:,2)), squeeze(TFR_spectra_sem(1,1,:,1)), 'b');
hold on;
shadedErrorBar(freqaxis, squeeze(TFR_spectra(2,1,:,2)), squeeze(TFR_spectra_sem(2,1,:,2)), 'r');
shadedErrorBar(freqaxis, squeeze(TFR_spectra(3,1,:,2)), squeeze(TFR_spectra_sem(3,1,:,2)), 'g');
shadedErrorBar(freqaxis, squeeze(TFR_spectra(4,1,:,2)), squeeze(TFR_spectra_sem(4,1,:,2)), 'm');
set(gca, 'xscale', 'log')
yline(0, 'k--')
xlim([freqaxis(1) freqaxis(end)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Perform cluster-based randomizations

%Combined monkeys
        low_current = [0.5 1.0];
        %MAry
    current_trials1 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 2);

%Jones
        low_current = [0.5 NaN];        
    current_trials2 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == low_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 1);
    
current_trials_low = [current_trials1; current_trials2];

%Combined monkeys
        high_current = [1.5 2.5];
%MAry
    current_trials1 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 2);

%Jones
        high_current = [1.0 1.5];
        
    current_trials2 = find((ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(1) | ALL_TFRbs_sessions{1}.trialinfo(:,2) == high_current(2)) & ...
        ALL_TFRbs_sessions{1}.trialinfo(:,5) == 1);
    
current_trials_high = [current_trials1; current_trials2];

% current_trials_high

spectensor = ALL_TFRbs_sessions;
spectensor{1}.powspctrm = spectensor{1}.powspctrm(current_trials_high,:,:,:);
spectensor{2}.powspctrm = spectensor{2}.powspctrm(current_trials_high,:,:,:);
spectensor{3}.powspctrm = spectensor{3}.powspctrm(current_trials_high,:,:,:);
spectensor{4}.powspctrm = spectensor{4}.powspctrm(current_trials_high,:,:,:);

% spectensor{1}.powspctrm = spectensor{1}.powspctrm(current_trials_low,:,:,:);
% spectensor{2}.powspctrm = spectensor{2}.powspctrm(current_trials_low,:,:,:);
% spectensor{3}.powspctrm = spectensor{3}.powspctrm(current_trials_low,:,:,:);
% spectensor{4}.powspctrm = spectensor{4}.powspctrm(current_trials_low,:,:,:);

rmindx = find(isnan(spectensor{1}.powspctrm(:,1,1,1)));
spectensor{1}.powspctrm(rmindx,:,:,:) = [];
spectensor{2}.powspctrm(rmindx,:,:,:) = [];
spectensor{3}.powspctrm(rmindx,:,:,:) = [];
spectensor{4}.powspctrm(rmindx,:,:,:) = [];


%Channel level randomization
arealist = {'PFC', '8A', 'PPC', 'STG'};
numfreq=length(freqaxis);
numtime = length(time_axis_baseline);
first_level_crit = 0.01;
corrected_alphaval = 0.01;

All_sigclust = zeros( numfreq, numtime,4);
All_empirical_effect = zeros(  numfreq, numtime,4);
nrand = 1000;
rand_clust_size = zeros(4,nrand);
for a = 1:4
    
    All_empirical_effect(:,:,a) = squeeze(mean(spectensor{a}.powspctrm(:,1,:,:)));
    
    baseline_bins = [nearest(time_axis_baseline, -10.0) nearest(time_axis_baseline, -3.0)];
    
    bs1 = nearest(time_axis_baseline,-10.0)
    bs2 = nearest(time_axis_baseline,-3);
    num_bs = length(bs1:bs2);
    during_DBS_bins = [nearest(time_axis_baseline, 3.0):nearest(time_axis_baseline, 28.5-3)];
    post_DBS_bins = [nearest(time_axis_baseline, 28.5+3):nearest(time_axis_baseline, 135.0)];
    num_post = length(during_DBS_bins)+length(post_DBS_bins);
    %perform randomization
    for r = 1:nrand
        TFR_rand = squeeze(spectensor{a}.powspctrm(:,:,:,:)); 
        numtrial = size(TFR_rand,1);
        for tr = 1:numtrial
            %randomize frequency and pre vs. post laser power
            baseline_bins = repmat(1:num_bs, [1 20]);
            baseline_bins = baseline_bins(randperm(length(baseline_bins)));
            baseline_bins = baseline_bins(1:num_post);
            drug_bins = [during_DBS_bins post_DBS_bins];
            allbins = [baseline_bins drug_bins];
            allbins = allbins(randperm(length(allbins)));
            rand_bins = allbins(1:numtime);

            TFR_rand(tr,:,:) = TFR_rand(tr,:,rand_bins);
        end
        
        %set artifact times to NaN
        artifact_bins = [nearest(time_axis_baseline, -3.0):nearest(time_axis_baseline, 3.0) nearest(time_axis_baseline, 28.5-3):nearest(time_axis_baseline, 28.5+3)];
        TFR_rand(:,:,artifact_bins) = NaN;
        first_level_stat = zeros(numfreq,numtime);
        baselineMAT = squeeze(nanmean(TFR_rand(:,:,bs1:bs2),3));
        baselineMAT = repmat(baselineMAT, [1 1 numtime]);
        
        [h p] = ttest2(baselineMAT, TFR_rand(:,:,:));
        first_level_stat(:,:) = squeeze(p);
        
        baselineMAT = squeeze(nanmean(TFR_rand(:,:,bs1:bs2),3));
        for t = 1:length(time_axis_baseline)
            [h p] = ttest2(baselineMAT, squeeze(TFR_rand(:,:,t)));
            first_level_stat(:,t) = squeeze(p);
        end

%         figure; imagesc(first_level_stat)
%         set(gca, 'ydir', 'normal')
%         colormap('Jet')
%         caxis([0 0.05])
% 
%         figure; imagesc(squeeze(mean(TFR_rand)))
%         set(gca, 'ydir', 'normal')
%         colormap('Jet')
%         caxis([-5 5])        
    
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
    figure; hist(rand_clust_size(a,:),100); title('cluster')
    
    sorted_clusters = sort(rand_clust_size(a,:), 'descend');
    alphabin = floor(corrected_alphaval.*nrand);
    threshold = sorted_clusters(alphabin);


    artifact_bins = [nearest(time_axis_baseline, -3.0):nearest(time_axis_baseline, 3.0) nearest(time_axis_baseline, 28.5-3):nearest(time_axis_baseline, 28.5+3)];
    TFR_empirical = squeeze(spectensor{a}.powspctrm(:,:,:,:));
    TFR_empirical(:,:,artifact_bins) = NaN;    
    
    first_level_stat = zeros(numfreq,numtime);
    %get empirical results
    for tbin = 1:numtime
        [h p] = ttest2(squeeze(nanmean(TFR_empirical(:,:,bs1:bs2),3)), TFR_empirical(:,:,tbin));
        first_level_stat(:,tbin) = p;
%         tbin
    end
    
    empirical_effect =  squeeze(first_level_stat(:,:)) < first_level_crit;
    empirical_mat_pass = zeros(numfreq, numtime);
    empirical_mat_pass(find(empirical_effect==1)) = 1;
    
    [L,num] = spm_bwlabel(empirical_mat_pass,6);
    Lpass = 0.3.*ones(size(L));
    Lpass_uncorrected = 0.3.*ones(size(L));
    tmp = find(first_level_stat<0.01);
    Lpass_uncorrected(tmp) = 1;
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
    All_empirical_effect(:,:,a) = squeeze(mean(TFR_empirical(:,:,:)));
    
    figure; imagesc(squeeze(mean(TFR_empirical(:,:,:))));
    set(gca, 'ydir', 'normal')
    caxis([-5 5])
    colorbar;
    alpha(Lpass_uncorrected)
    hold on;
    xline(nearest(time_axis_baseline, 0), 'k--')
    title(arealist{a})
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    xlim([xlim1 xlim2])
    xlabel('Time from DBS onset (0 s)');
    ylabel('Frequency (Hz)')    
    colormap('Jet')
    title(['Uncorrected - ' arealist{a}])
    
    figure; imagesc(squeeze(mean(TFR_empirical(:,:,:))));
    set(gca, 'ydir', 'normal')
    caxis([-5 5])
    colorbar;
    alpha(Lpass)
    hold on;
    xline(nearest(time_axis_baseline, 0), 'k--')
    title(arealist{a})
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    xlim([xlim1 xlim2])
    xlabel('Time from DBS onset (0 s)');
    ylabel('Frequency (Hz)')    
    colormap('Jet')
    title(['Corrected - ' arealist{a}])
    
    pause(0.1)
    
    
    
end
    



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Low current

spectensor = ALL_TFRbs_sessions;
% spectensor{1}.powspctrm = spectensor{1}.powspctrm(current_trials_high,:,:,:);
% spectensor{2}.powspctrm = spectensor{2}.powspctrm(current_trials_high,:,:,:);
% spectensor{3}.powspctrm = spectensor{3}.powspctrm(current_trials_high,:,:,:);
% spectensor{4}.powspctrm = spectensor{4}.powspctrm(current_trials_high,:,:,:);

spectensor{1}.powspctrm = spectensor{1}.powspctrm(current_trials_low,:,:,:);
spectensor{2}.powspctrm = spectensor{2}.powspctrm(current_trials_low,:,:,:);
spectensor{3}.powspctrm = spectensor{3}.powspctrm(current_trials_low,:,:,:);
spectensor{4}.powspctrm = spectensor{4}.powspctrm(current_trials_low,:,:,:);

rmindx = find(isnan(spectensor{1}.powspctrm(:,1,1,1)));
spectensor{1}.powspctrm(rmindx,:,:,:) = [];
spectensor{2}.powspctrm(rmindx,:,:,:) = [];
spectensor{3}.powspctrm(rmindx,:,:,:) = [];
spectensor{4}.powspctrm(rmindx,:,:,:) = [];


%Channel level randomization
arealist = { 'PFC', '8A' 'PPC', 'STG'};
numfreq=length(freqaxis);
numtime = length(time_axis_baseline);
first_level_crit = 0.01;
corrected_alphaval = 0.01;

All_sigclust = zeros( numfreq, numtime,4);
All_empirical_effect = zeros(  numfreq, numtime,4);
nrand = 1000;
rand_clust_size = zeros(4,nrand);
for a = 1:4
    
    All_empirical_effect(:,:,a) = squeeze(mean(spectensor{a}.powspctrm(:,1,:,:)));
    
    baseline_bins = [nearest(time_axis_baseline, -10.0) nearest(time_axis_baseline, -3.0)];
    
    bs1 = nearest(time_axis_baseline,-10.0)
    bs2 = nearest(time_axis_baseline,-3);
    num_bs = length(bs1:bs2);
    during_DBS_bins = [nearest(time_axis_baseline, 3.0):nearest(time_axis_baseline, 28.5-3)];
    post_DBS_bins = [nearest(time_axis_baseline, 28.5+3):nearest(time_axis_baseline, 135.0)];
    num_post = length(during_DBS_bins)+length(post_DBS_bins);
    %perform randomization
    for r = 1:nrand
        TFR_rand = squeeze(spectensor{a}.powspctrm(:,:,:,:)); 
        numtrial = size(TFR_rand,1);
        for tr = 1:numtrial
            %randomize frequency and pre vs. post laser power
            baseline_bins = repmat(1:num_bs, [1 20]);
            baseline_bins = baseline_bins(randperm(length(baseline_bins)));
            baseline_bins = baseline_bins(1:num_post);
            drug_bins = [during_DBS_bins post_DBS_bins];
            allbins = [baseline_bins drug_bins];
            allbins = allbins(randperm(length(allbins)));
            rand_bins = allbins(1:numtime);

            TFR_rand(tr,:,:) = TFR_rand(tr,:,rand_bins);
        end
        
        %set artifact times to NaN
        artifact_bins = [nearest(time_axis_baseline, -3.0):nearest(time_axis_baseline, 3.0) nearest(time_axis_baseline, 28.5-3):nearest(time_axis_baseline, 28.5+3)];
        TFR_rand(:,:,artifact_bins) = NaN;
        first_level_stat = zeros(numfreq,numtime);
        baselineMAT = squeeze(nanmean(TFR_rand(:,:,bs1:bs2),3));
        baselineMAT = repmat(baselineMAT, [1 1 numtime]);
        
        [h p] = ttest2(baselineMAT, TFR_rand(:,:,:));
        first_level_stat(:,:) = squeeze(p);
        
        baselineMAT = squeeze(nanmean(TFR_rand(:,:,bs1:bs2),3));
        for t = 1:length(time_axis_baseline)
            [h p] = ttest2(baselineMAT, squeeze(TFR_rand(:,:,t)));
            first_level_stat(:,t) = squeeze(p);
        end

%         figure; imagesc(first_level_stat)
%         set(gca, 'ydir', 'normal')
%         colormap('Jet')
%         caxis([0 0.05])
% 
%         figure; imagesc(squeeze(mean(TFR_rand)))
%         set(gca, 'ydir', 'normal')
%         colormap('Jet')
%         caxis([-5 5])        
    
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
    figure; hist(rand_clust_size(a,:),100); title('cluster')
    
    sorted_clusters = sort(rand_clust_size(a,:), 'descend');
    alphabin = floor(corrected_alphaval.*nrand);
    threshold = sorted_clusters(alphabin);


    artifact_bins = [nearest(time_axis_baseline, -3.0):nearest(time_axis_baseline, 3.0) nearest(time_axis_baseline, 28.5-3):nearest(time_axis_baseline, 28.5+3)];
    TFR_empirical = squeeze(spectensor{a}.powspctrm(:,:,:,:));
    TFR_empirical(:,:,artifact_bins) = NaN;    
    
    first_level_stat = zeros(numfreq,numtime);
    %get empirical results
    for tbin = 1:numtime
        [h p] = ttest2(squeeze(nanmean(TFR_empirical(:,:,bs1:bs2),3)), TFR_empirical(:,:,tbin));
        first_level_stat(:,tbin) = p;
%         tbin
    end
    
    empirical_effect =  squeeze(first_level_stat(:,:)) < first_level_crit;
    empirical_mat_pass = zeros(numfreq, numtime);
    empirical_mat_pass(find(empirical_effect==1)) = 1;
    
    [L,num] = spm_bwlabel(empirical_mat_pass,6);
    Lpass = 0.3.*ones(size(L));
    Lpass_uncorrected = 0.3.*ones(size(L));
    tmp = find(first_level_stat<0.01);
    Lpass_uncorrected(tmp) = 1;
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
    All_empirical_effect(:,:,a) = squeeze(mean(TFR_empirical(:,:,:)));
    
    figure; imagesc(squeeze(mean(TFR_empirical(:,:,:))));
    set(gca, 'ydir', 'normal')
    caxis([-5 5])
    colorbar;
    alpha(Lpass_uncorrected)
    hold on;
    xline(nearest(time_axis_baseline, 0), 'k--')
    title(arealist{a})
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    xlim([xlim1 xlim2])
    xlabel('Time from DBS onset (0 s)');
    ylabel('Frequency (Hz)')    
    colormap('Jet')
    title(['Uncorrected - ' arealist{a}])
    
    figure; imagesc(squeeze(mean(TFR_empirical(:,:,:))));
    set(gca, 'ydir', 'normal')
    caxis([-5 5])
    colorbar;
    alpha(Lpass)
    hold on;
    xline(nearest(time_axis_baseline, 0), 'k--')
    title(arealist{a})
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
    xlim([xlim1 xlim2])
    xlabel('Time from DBS onset (0 s)');
    ylabel('Frequency (Hz)')    
    colormap('Jet')
    title(['Corrected - ' arealist{a}])
    
    pause(0.1)
    
    
    
end
    

cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')

save DBS_TFR_stats_1000rand_low_current All_sigclust All_empirical_effect xmarkers xlabels ymarkers ylabels arealist time_axis_baseline





%%
%load up data from file
cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')
load DBS_TFR_stats_1000rand_high_current 
arealist = { 'PFC', '8A' 'PPC', 'STG'};

    %create time axis
    xmarkers = [nearest(time_axis_baseline, -30) nearest(time_axis_baseline, 0) nearest(time_axis_baseline, 30) ...
        nearest(time_axis_baseline, 60) nearest(time_axis_baseline, 90), nearest(time_axis_baseline, 120)];
    
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis_baseline(xmarkers(indx)));
        indx=indx+1;
    end
    
a=4;
    figure; imagesc(squeeze(All_empirical_effect(:,:,a)));
    set(gca, 'ydir', 'normal')
    caxis([-8 8])
    colorbar;
    alpha(squeeze(All_sigclust(:,:,a)))
    hold on;
    contourLines = zeros(size(squeeze(All_sigclust(:,:,a))));
    contourLines(find(squeeze(All_sigclust(:,:,a))==1)) = 1;
    contour(contourLines,1, 'k')    
    xline(nearest(time_axis_baseline, 0), 'k--')
    title(arealist{a})
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
%     xlim([xlim1 xlim2])
    xlabel('Time from DBS onset (0 s)');
    ylabel('Frequency (Hz)')    
    colormap('Jet')
    title(['High Current - ' arealist{a}])


load DBS_TFR_stats_1000rand_low_current 
arealist = { 'PFC', '8A' 'PPC', 'STG'};


    %create time axis
    xmarkers = [nearest(time_axis_baseline, -30) nearest(time_axis_baseline, 0) nearest(time_axis_baseline, 30) ...
        nearest(time_axis_baseline, 60) nearest(time_axis_baseline, 90), nearest(time_axis_baseline, 120)];
    
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(time_axis_baseline(xmarkers(indx)));
        indx=indx+1;
    end

    figure; imagesc(squeeze(All_empirical_effect(:,:,a)));
    set(gca, 'ydir', 'normal')
    caxis([-8 8])
    colorbar;
    alpha(squeeze(All_sigclust(:,:,a)))
    hold on;
    contourLines = zeros(size(squeeze(All_sigclust(:,:,a))));
    contourLines(find(squeeze(All_sigclust(:,:,a))==1)) = 1;
    contour(contourLines,1, 'k')  
    xline(nearest(time_axis_baseline, 0), 'k--')
    title(arealist{a})
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    set(gca, 'ytick', ymarkers)
    set(gca, 'yticklabel', ylabels)
%     xlim([xlim1 xlim2])
    xlabel('Time from DBS onset (0 s)');
    ylabel('Frequency (Hz)')    
    colormap('Jet')
    title(['Low Current - ' arealist{a}])




















    
    
    
    