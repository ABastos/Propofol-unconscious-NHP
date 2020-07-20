
session_list = ls('\\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolPuffTone\*.mat');
session_list(end,:) = []

for s = 2:size(session_list,1)
    cd('\\millerdata.mit.edu\common\datasets\anesthesia\mat\propofolPuffTone\')
    
    load(deblank(session_list(s,:)), 'spikeTimes', 'unitInfo', 'sessionInfo', 'trialInfo')
    
    %events in seconds
%     ntime = size(lfp,1);
    timevec_seconds= -(20*60):5:(10*60);
    ROC = sessionInfo.eyesOpen(1);
    LOC = sessionInfo.eyesClose(1)-ROC;
    drugStart = sessionInfo.drugStart(1)-ROC;
    drugStart2 = sessionInfo.drugStart(2)-ROC;
    drugEnd = sessionInfo.drugEnd(2)-ROC;


    PFCindx = find(strcmp(unitInfo.array, 'PFC') & unitInfo.unit > 0);
    FEFindx = find(strcmp(unitInfo.array, 'FEF') & unitInfo.unit > 0);
    PPCindx = find(strcmp(unitInfo.array, 'PPC') & unitInfo.unit > 0);
    STGindx = find(strcmp(unitInfo.array, 'STG') & unitInfo.unit > 0);
    

    %get all stimulus times
    time_exclude = [];
    for tr = 1:length(trialInfo.cpt_puffOn)
        if ~isnan(trialInfo.cpt_puffOn(tr)) && ~isnan(trialInfo.cpt_toneOn(tr))
            time_exclude = [time_exclude; trialInfo.cpt_toneOn(tr) trialInfo.cpt_puffOff(tr)+0.25];
        elseif ~isnan(trialInfo.cpt_puffOn(tr)) && isnan(trialInfo.cpt_toneOn(tr))
            time_exclude = [time_exclude; trialInfo.cpt_puffOn(tr) trialInfo.cpt_puffOff(tr)+0.25];
        elseif ~isnan(trialInfo.cpt_toneOn(tr))  && isnan(trialInfo.cpt_puffOn(tr))
            time_exclude = [time_exclude; trialInfo.cpt_toneOn(tr) trialInfo.cpt_toneOff(tr)+0.25];
        end        
    end

    time_exclude = time_exclude - LOC;
    


    spike_rates = cell(4,1);
    spike_rates{1} = zeros(length(PFCindx),length(timevec_seconds));
    spike_rates{2} = zeros(length(FEFindx),length(timevec_seconds));
    spike_rates{3} = zeros(length(PPCindx),length(timevec_seconds));
    spike_rates{4} = zeros(length(STGindx),length(timevec_seconds));
    
    integration_window = 5;
    
    %remove times around puff/tone stimulation
    
    for t = 1:length(timevec_seconds)

            latency = [timevec_seconds(t)-integration_window/2:0.001:timevec_seconds(t)+integration_window/2];

            rmi = [];
            newtime = [];
            for exclude_tr = 1:size(time_exclude,1)
%                 if (latency(1) > time_exclude(exclude_tr,1) && latency(1) < time_exclude(exclude_tr,2)) || ...
%                         (latency(end) > time_exclude(exclude_tr,1) && latency(end) < time_exclude(exclude_tr,2))
                    
                rmi = [rmi intersect(latency(1):0.001:latency(end), time_exclude(exclude_tr,1):0.001:time_exclude(exclude_tr,2))];

%                     rmi = [rmi; exclude_tr];
%                     
%                     %find time to exclude
%                     time_exclude_tr = time_exclude(exclude_tr,1):0.001:time_exclude(exclude_tr,2);
%                     newtime = [newtime setdiff(latency, time_exclude_tr)];
%                     
%                 end
            end
            
            latency = setdiff(latency, rmi);
            sniplen = length(latency)/1000;
            %PFC
            for c = 1:length(PFCindx)
                spike_rates{1}(c,t) = length(find(spikeTimes{PFCindx(c)}-ROC<=latency(end) & spikeTimes{PFCindx(c)}-ROC>=latency(1)))/sniplen;
            end

            %FEF
            for c = 1:length(FEFindx)
                spike_rates{2}(c,t) = length(find(spikeTimes{FEFindx(c)}-ROC<=latency(end) & spikeTimes{FEFindx(c)}-ROC>=latency(1)))/sniplen;
            end
            %PPC
            for c = 1:length(PPCindx)
                spike_rates{3}(c,t) = length(find(spikeTimes{PPCindx(c)}-ROC<=latency(end) & spikeTimes{PPCindx(c)}-ROC>=latency(1)))/sniplen;
            end
            %STG
            for c = 1:length(STGindx)
                spike_rates{4}(c,t) = length(find(spikeTimes{STGindx(c)}-ROC<=latency(end) & spikeTimes{STGindx(c)}-ROC>=latency(1)))/sniplen;
            end
            t
            
    end        

    %create time axis
    xmarkers = 18:25:length(timevec_seconds);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(timevec_seconds(xmarkers(indx))/60);
        indx=indx+1;
    end 
    
    figure('Name', sessionInfo.session, 'units','normalized','outerposition',[0 0 1 1]) 
    subplot(2,2,1);
    imagesc(spike_rates{1}); hold on; colorbar
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    ylabel('Unit index')
    xline(nearest(timevec_seconds, 0), 'k', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, LOC), 'k', 'LineWidth', 2)
    xline(nearest(timevec_seconds, drugEnd), 'w', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, ROC), 'w', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, drugStart), 'w', 'LineWidth', 2)
    title('PFC firing rates')
    caxis([0 20])
    
    subplot(2,2,2);
    imagesc(spike_rates{2}); hold on; colorbar
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    ylabel('Unit index')
    xline(nearest(timevec_seconds, 0), 'k', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, LOC), 'k', 'LineWidth', 2)
    xline(nearest(timevec_seconds, drugEnd), 'w', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, ROC), 'w', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, drugStart), 'w', 'LineWidth', 2)
    title('FEF firing rates')
    caxis([0 10])
    
    subplot(2,2,3); 
    imagesc(spike_rates{3}); hold on; colorbar
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    ylabel('Unit index')
    xline(nearest(timevec_seconds, 0), 'k', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, LOC), 'k', 'LineWidth', 2)
    xline(nearest(timevec_seconds, drugEnd), 'w', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, ROC), 'w', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, drugStart), 'w', 'LineWidth', 2)
    title('PPC firing rates')
    caxis([0 10])
    
    subplot(2,2,4); 
    imagesc(spike_rates{4}); hold on; colorbar
    set(gca, 'ydir', 'normal')
    colormap('Jet')
    set(gca, 'xtick', xmarkers)
    set(gca, 'xticklabel', xlabels)
    xlabel('Time since Drug start (minutes)')
    ylabel('Unit index')
    xline(nearest(timevec_seconds, 0), 'k', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, LOC), 'k', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, drugStart), 'w', 'LineWidth', 2)
    xline(nearest(timevec_seconds, drugEnd), 'w', 'LineWidth', 2)
%     xline(nearest(timevec_seconds, ROC), 'w', 'LineWidth', 2)
    title('STG firing rates')
    caxis([0 10])
    
%     figure; shadedErrorBar(timevec_seconds/60, mean(spike_rates{1}), std(spike_rates{1})./sqrt(length(spike_rates{1})))
%     xlim([-15 75])
%     xline(0, 'k', 'LineWidth', 2)
%     xline(LOC, 'k', 'LineWidth', 2)
%     xline(drugEnd, 'r', 'LineWidth', 2)
%     xline(ROC, 'g', 'LineWidth', 2)
%     xlabel('Time since Drug start (minutes)')
%     ylabel('Firing rates')
%     title('PFC')
%     
%     figure; shadedErrorBar(timevec_seconds/60, mean(spike_rates{2}), std(spike_rates{2})./sqrt(length(spike_rates{2})))
%     xlim([-15 75])
%     xline(0, 'k', 'LineWidth', 2)
%     xline(LOC, 'k', 'LineWidth', 2)
%     xline(drugEnd, 'r', 'LineWidth', 2)
%     xline(ROC, 'g', 'LineWidth', 2)
%     xlabel('Time since Drug start (minutes)')
%     ylabel('Firing rates')
%     title('FEF')
%     
%     figure; shadedErrorBar(timevec_seconds/60, mean(spike_rates{3}), std(spike_rates{3})./sqrt(length(spike_rates{3})))
%     xlim([-15 75])
%     xline(0, 'k', 'LineWidth', 2)
%     xline(LOC, 'k', 'LineWidth', 2)
%     xline(drugEnd, 'r', 'LineWidth', 2)
%     xline(ROC, 'g', 'LineWidth', 2)
%     xlabel('Time since Drug start (minutes)')
%     ylabel('Firing rates')
%     title('PPC')
%     
%     figure; shadedErrorBar(timevec_seconds/60, mean(spike_rates{4}), std(spike_rates{4})./sqrt(length(spike_rates{4})))
%     xlim([-15 75])
%     xline(0, 'k', 'LineWidth', 2)
%     xline(LOC, 'k', 'LineWidth', 2)
%     xline(drugEnd, 'r', 'LineWidth', 2)
%     xline(ROC, 'g', 'LineWidth', 2)
%     xlabel('Time since Drug start (minutes)')
%     ylabel('Firing rates')
%     title('STG')
    
    
    figure('Name', sessionInfo.session); 
    shadedErrorBar(timevec_seconds/60, mean(spike_rates{1},1), std(spike_rates{1},[],1)./sqrt(length(spike_rates{1})), 'g'); hold on
    shadedErrorBar(timevec_seconds/60, mean(spike_rates{2},1), std(spike_rates{2},[],1)./sqrt(length(spike_rates{2})), 'r')
    shadedErrorBar(timevec_seconds/60, mean(spike_rates{3},1), std(spike_rates{3},[],1)./sqrt(length(spike_rates{3})), 'm')
    shadedErrorBar(timevec_seconds/60, mean(spike_rates{4},1), std(spike_rates{4},[],1)./sqrt(length(spike_rates{4})), 'b')
    xlim([-30 10])
    xline(0, 'k', 'LineWidth', 2)
%     xline(LOC, 'k', 'LineWidth', 2)
    xline(drugEnd/60, 'r', 'LineWidth', 2)
%     xline(ROC, 'g', 'LineWidth', 2)
%     xline(drugStart/60, 'g', 'LineWidth', 2)
    xlabel('Time since Drug start (minutes)')
    ylabel('Firing rates')
    
    
    
    cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')
    save([sessionInfo.session '_SUA_ROC_Locked_NoPuffTones'], 'drugStart', 'drugStart2', 'drugEnd', 'LOC', 'ROC', 'timevec_seconds', 'sessionInfo', 'spike_rates', '-v7.3')
    
    clear data lfp TFR TFRbs sessionInfo spike_rates
    
end

return

%load in data for average spectrogram and spiking

%%

cd('\\millerdata.mit.edu\common\Andre\Anesthesia_Paper_Analysis\')


filestoload = ls('*_SUA_ROC_Locked_NoPuffTones*')
timevec_seconds= -(20*60):5:(10*60);


spike_rates_sessions = {};

nunits_global = ones(4,1);
drugstart_sessions = zeros(size(filestoload,1),1);
drugstart2_sessions = zeros(size(filestoload,1),1);
drugstop_sessions = zeros(size(filestoload,1),1);
LOC_sessions = zeros(size(filestoload,1),1);
ROC_sessions = zeros(size(filestoload,1),1);
for f = 1:size(filestoload,1)

    load(filestoload(f,:)); %, 'sessionInfo', 'spike_rates', 'timevec_seconds')

    
    for a = 1:4
        nunits = size(spike_rates{a},1);
        if f ==1
            spike_rates_sessions{a} = spike_rates{a};
        else
            spike_rates_sessions{a} = cat(1, spike_rates_sessions{a}, spike_rates{a});
        end
        nunits_global(a) = nunits_global(a) + nunits;
    end
    
%         
    drugstart_sessions(f) = drugStart/60;
    drugstart2_sessions(f) = drugStart2/60;
    drugstop_sessions(f) = drugEnd/60;
    LOC_sessions(f) = LOC/60;
    ROC_sessions(f) = ROC/60;
    
       
    f
    
end

%%

timevec_minutes = timevec_seconds/60;
    %create time axis
    xmarkers = 18:25:length(timevec_minutes);
    xlabels = {};
    indx=1;
    for x = xmarkers
        xlabels{indx} = int2str(timevec_minutes(xmarkers(indx)));
        indx=indx+1;
    end 
    


spike_rates_CAT_ROC = cat(1, spike_rates_sessions{1}, spike_rates_sessions{2}, spike_rates_sessions{3},spike_rates_sessions{4});

ROCbin = nearest(timevec_minutes, 0)
meanFRROC = mean(mean(spike_rates_CAT_ROC(:,ROCbin),2));
semFRROC = std(mean(spike_rates_CAT_ROC(:,ROCbin),2))./size(spike_rates_CAT_ROC,1);


mean(mean(spike_rates_sessions{1}(:,ROCbin),2))
mean(mean(spike_rates_sessions{2}(:,ROCbin),2))
mean(mean(spike_rates_sessions{3}(:,ROCbin),2))
mean(mean(spike_rates_sessions{4}(:,ROCbin),2))


ROCbin1 = nearest(timevec_minutes, -20);
ROCbin2 = nearest(timevec_minutes, -15);
meanFRROC = mean(mean(spike_rates_CAT_ROC(:,ROCbin1:ROCbin2),2));
semFRROC = std(mean(spike_rates_CAT_ROC(:,ROCbin1:ROCbin2),2))./size(spike_rates_CAT_ROC,1);



mean(mean(spike_rates_sessions{1}(:,ROCbin1:ROCbin2),2))
mean(mean(spike_rates_sessions{2}(:,ROCbin1:ROCbin2),2))
mean(mean(spike_rates_sessions{3}(:,ROCbin1:ROCbin2),2))
mean(mean(spike_rates_sessions{4}(:,ROCbin1:ROCbin2),2))

spike_sem = zeros(4, length(timevec_seconds));
spike_sem(1,:) = nanstd(spike_rates_sessions{1})./sqrt(size(spike_rates_sessions{1},1));
spike_sem(2,:) = nanstd(spike_rates_sessions{2})./sqrt(size(spike_rates_sessions{2},1));
spike_sem(3,:) = nanstd(spike_rates_sessions{3})./sqrt(size(spike_rates_sessions{3},1));
spike_sem(4,:) = nanstd(spike_rates_sessions{4})./sqrt(size(spike_rates_sessions{4},1));

%downsample spiking
newtimevec =  downsample(timevec_minutes,6); %one measurement every 30 seconds
spikes_smdn = {};
spikes_smdn{1} = zeros(size(spike_rates_sessions{1},1), length(newtimevec));
spikes_smdn{2} = zeros(size(spike_rates_sessions{2},1), length(newtimevec));
spikes_smdn{3} = zeros(size(spike_rates_sessions{3},1), length(newtimevec));
spikes_smdn{4} = zeros(size(spike_rates_sessions{4},1), length(newtimevec));
for a = 1:4
    for c = 1:size(spike_rates_sessions{a},1)
        spikes_smdn{a}(c,:) = downsample(smooth(spike_rates_sessions{a}(c,:),6),6);
    end
    a
end

%calculate the bootstrap CI for each time bin
spike_CI = zeros(4, 2, length(newtimevec)); %upper/lower 99% CI
for a = 1:4
    for t = 1:length(newtimevec)
        [confInts] = bootstrapConfIntervals(spikes_smdn{a}(:,t), @nanmean, [],0.99,10000);
        spike_CI(a,1,t) = confInts(1) - nanmean(spikes_smdn{a}(:,t)); %upper bound
        spike_CI(a,2,t) = nanmean(spikes_smdn{a}(:,t)) - confInts(2); %lower bound
        t
    end
end

figure; 
shadedErrorBar(newtimevec, nanmean(spikes_smdn{1}), squeeze(spike_CI(1,:,:)), 'b'); hold on;
shadedErrorBar(newtimevec, nanmean(spikes_smdn{2}), squeeze(spike_CI(2,:,:)), 'r'); hold on;
shadedErrorBar(newtimevec, nanmean(spikes_smdn{3}), squeeze(spike_CI(3,:,:)), 'g'); hold on;
shadedErrorBar(newtimevec, nanmean(spikes_smdn{4}), squeeze(spike_CI(4,:,:)), 'm'); hold on;
xline(0, 'k--', 'LineWidth', 1)
xline(mean(drugstop_sessions), 'k--', 'LineWidth', 1)
xline(mean(LOC_sessions), 'r', 'LineWidth', 2)
xline(mean(LOC_sessions)-std(LOC_sessions)./sqrt(size(LOC_sessions,1)), 'r', 'LineWidth', 1)
xline(mean(LOC_sessions)+std(LOC_sessions)./sqrt(size(LOC_sessions,1)), 'r', 'LineWidth', 1)
xline(mean(drugstop_sessions), 'k', 'LineWidth', 2)
xline(mean(drugstop_sessions)-std(drugstop_sessions)./sqrt(size(drugstop_sessions,1)), 'k', 'LineWidth', 1)
xline(mean(drugstop_sessions)+std(drugstop_sessions)./sqrt(size(drugstop_sessions,1)), 'k', 'LineWidth', 1)
xlim([-20 10])
xlabel('Time since ROC (minutes)')
ylabel('Firing rates (spikes per second)')
title('PFC=blue, 8A=red, PPC=green, STG=magenta')


ROCbin = nearest(newtimevec, 0)
meanFRROC = zeros(4,1);
meanFRROC(1) = nanmean(spikes_smdn{1}(:,ROCbin))
meanFRROC(2) = nanmean(spikes_smdn{2}(:,ROCbin))
meanFRROC(3) = nanmean(spikes_smdn{3}(:,ROCbin))
meanFRROC(4) = nanmean(spikes_smdn{4}(:,ROCbin))

meanCIFRLOC = zeros(4,2);
[confInts] = bootstrapConfIntervals(spikes_smdn{1}(:,ROCbin), @nanmean, [],0.99,10000);
meanCIFRLOC(1,1) = confInts(1); 
meanCIFRLOC(1,2) = confInts(2); 
[confInts] = bootstrapConfIntervals(spikes_smdn{2}(:,ROCbin), @nanmean, [],0.99,10000);
meanCIFRLOC(2,1) = confInts(1); 
meanCIFRLOC(2,2) = confInts(2); 
[confInts] = bootstrapConfIntervals(spikes_smdn{3}(:,ROCbin), @nanmean, [],0.99,10000);
meanCIFRLOC(3,1) = confInts(1); 
meanCIFRLOC(3,2) = confInts(2); 
[confInts] = bootstrapConfIntervals(spikes_smdn{4}(:,ROCbin), @nanmean, [],0.99,10000);
meanCIFRLOC(4,1) = confInts(1); 
meanCIFRLOC(4,2) = confInts(2); 

figure; 
shadedErrorBar(downsample(timevec_minutes,6), downsample(mean(spike_rates_sessions{4}),6), downsample(spike_sem(4,:),6), 'b'); hold on;
shadedErrorBar(timevec_minutes, mean(spike_rates_sessions{2}), spike_sem(2,:), 'r'); hold on;
shadedErrorBar(timevec_minutes, mean(spike_rates_sessions{3}), spike_sem(3,:), 'g'); hold on;
shadedErrorBar(timevec_minutes, mean(spike_rates_sessions{4}), spike_sem(4,:), 'm'); hold on;
xline(0, 'k--', 'LineWidth', 1)
xline(mean(drugstop_sessions), 'k--', 'LineWidth', 1)
xline(mean(LOC_sessions), 'r', 'LineWidth', 2)
xline(mean(LOC_sessions)-std(LOC_sessions)./sqrt(size(LOC_sessions,1)), 'r', 'LineWidth', 1)
xline(mean(LOC_sessions)+std(LOC_sessions)./sqrt(size(LOC_sessions,1)), 'r', 'LineWidth', 1)
% xline(mean(ROC_sessions), 'g', 'LineWidth', 2)
% xline(mean(ROC_sessions)-std(ROC_sessions)./sqrt(size(ROC_sessions,1)), 'g', 'LineWidth', 1)
% xline(mean(ROC_sessions)+std(ROC_sessions)./sqrt(size(ROC_sessions,1)), 'g', 'LineWidth', 1)
xline(mean(drugstop_sessions), 'k', 'LineWidth', 2)
xline(mean(drugstop_sessions)-std(drugstop_sessions)./sqrt(size(drugstop_sessions,1)), 'k', 'LineWidth', 1)
xline(mean(drugstop_sessions)+std(drugstop_sessions)./sqrt(size(drugstop_sessions,1)), 'k', 'LineWidth', 1)
% xline(mean(drugstart2_sessions), 'g', 'LineWidth', 2)
% xline(mean(drugstart2_sessions)-std(drugstart2_sessions)./sqrt(size(drugstart2_sessions,1)), 'g', 'LineWidth', 1)
% xline(mean(drugstart2_sessions)+std(drugstart2_sessions)./sqrt(size(drugstart2_sessions,1)), 'g', 'LineWidth', 1)
xlim([-20 10])
xlabel('Time since Drug start (minutes)')
ylabel('Firing rates (spikes per second)')
title('PFC=blue, 8A=red, PPC=green, STG=magenta')
% set(gca, 'yscale', 'log')
%%

ROCbin = nearest(timevec_minutes, 0)

mean(mean(spike_rates_sessions{1}(:,ROCbin),2))
mean(mean(spike_rates_sessions{2}(:,ROCbin),2))
mean(mean(spike_rates_sessions{3}(:,ROCbin),2))
mean(mean(spike_rates_sessions{4}(:,ROCbin),2))

[p h ] = ranksum( mean(spike_rates_sessions{1}(:,ROCbin),2), mean(spike_rates_sessions{3}(:,ROCbin),2))
[p h ] = ranksum( mean(spike_rates_sessions{1}(:,ROCbin),2), mean(spike_rates_sessions{4}(:,ROCbin),2))
[p h ] = ranksum( mean(spike_rates_sessions{2}(:,ROCbin),2), mean(spike_rates_sessions{3}(:,ROCbin),2))
[p h ] = ranksum( mean(spike_rates_sessions{2}(:,ROCbin),2), mean(spike_rates_sessions{4}(:,ROCbin),2))


ROCbin = nearest(timevec_minutes, mean(ROC_sessions))

mean(mean(spike_rates_sessions{1}(:,ROCbin),2))
mean(mean(spike_rates_sessions{2}(:,ROCbin),2))
mean(mean(spike_rates_sessions{3}(:,ROCbin),2))
mean(mean(spike_rates_sessions{4}(:,ROCbin),2))

[p h ] = ranksum( mean(spike_rates_sessions{1}(:,ROCbin),2), mean(spike_rates_sessions{3}(:,ROCbin),2))
[p h ] = ranksum( mean(spike_rates_sessions{1}(:,ROCbin),2), mean(spike_rates_sessions{4}(:,ROCbin),2))
[p h ] = ranksum( mean(spike_rates_sessions{2}(:,ROCbin),2), mean(spike_rates_sessions{3}(:,ROCbin),2))
[p h ] = ranksum( mean(spike_rates_sessions{2}(:,ROCbin),2), mean(spike_rates_sessions{4}(:,ROCbin),2))

[p h ] = ranksum( mean(spike_rates_sessions{3}(:,ROCbin),2), mean(spike_rates_sessions{4}(:,ROCbin),2))

%%


bs1=nearest(timevec_minutes,-15);
bs2=nearest(timevec_minutes,0);

LOC1=nearest(timevec_minutes,30);
LOC2=nearest(timevec_minutes,60);

mean(mean(spike_rates_sessions{1}(:,bs1:bs2),2))

FR_awake = zeros(4,2);
FR_LOC = zeros(4,2);

FR_awake(1,1) = mean(mean(spike_rates_sessions{1}(:,bs1:bs2),2));
FR_awake(2,1) = mean(mean(spike_rates_sessions{2}(:,bs1:bs2),2));
FR_awake(3,1) = mean(mean(spike_rates_sessions{3}(:,bs1:bs2),2));
FR_awake(4,1) = mean(mean(spike_rates_sessions{4}(:,bs1:bs2),2));

FR_awake(1,2) = std(mean(spike_rates_sessions{1}(:,bs1:bs2),2))./sqrt(size(spike_rates_sessions{1},1));
FR_awake(2,2) = std(mean(spike_rates_sessions{2}(:,bs1:bs2),2))./sqrt(size(spike_rates_sessions{2},1));
FR_awake(3,2) = std(mean(spike_rates_sessions{3}(:,bs1:bs2),2))./sqrt(size(spike_rates_sessions{3},1));
FR_awake(4,2) = std(mean(spike_rates_sessions{4}(:,bs1:bs2),2))./sqrt(size(spike_rates_sessions{4},1));

FR_LOC(1,1) = mean(mean(spike_rates_sessions{1}(:,LOC1:LOC2),2));
FR_LOC(2,1) = mean(mean(spike_rates_sessions{2}(:,LOC1:LOC2),2));
FR_LOC(3,1) = mean(mean(spike_rates_sessions{3}(:,LOC1:LOC2),2));
FR_LOC(4,1) = mean(mean(spike_rates_sessions{4}(:,LOC1:LOC2),2));

FR_LOC(1,2) = std(mean(spike_rates_sessions{1}(:,LOC1:LOC2),2))./sqrt(size(spike_rates_sessions{1},1));
FR_LOC(2,2) = std(mean(spike_rates_sessions{2}(:,LOC1:LOC2),2))./sqrt(size(spike_rates_sessions{2},1));
FR_LOC(3,2) = std(mean(spike_rates_sessions{3}(:,LOC1:LOC2),2))./sqrt(size(spike_rates_sessions{3},1));
FR_LOC(4,2) = std(mean(spike_rates_sessions{4}(:,LOC1:LOC2),2))./sqrt(size(spike_rates_sessions{4},1));

for a=1:4
tmp = mean(spike_rates_sessions{a}(:,bs1:bs2),2);
rmi = find(tmp==0); 
spike_rates_sessions{a}(rmi,:) = [];
tmp = mean(spike_rates_sessions{a}(:,LOC1:LOC2),2);
rmi = find(tmp==0); 
spike_rates_sessions{a}(rmi,:) = [];
end

effect_size = zeros(4,1);
effect_size_sem = zeros(4,1);
for a =1:4
effect_size(a) = mean(10.*log10(mean(spike_rates_sessions{a}(:,LOC1:LOC2),2)./mean(spike_rates_sessions{a}(:,bs1:bs2),2)));
nunits = size(spike_rates_sessions{a},1);
effect_size_sem(a) = std(10.*log10(mean(spike_rates_sessions{a}(:,LOC1:LOC2),2)./mean(spike_rates_sessions{a}(:,bs1:bs2),2)))./sqrt(nunits);
end

%effect size comparison
effect_size_compare = NaN(4,4);
for a1=1:4 
    for a2=a1+1:4
        tmp1 = 10.*log10(mean(spike_rates_sessions{a1}(:,LOC1:LOC2),2)./mean(spike_rates_sessions{a1}(:,bs1:bs2),2));
        tmp2 = 10.*log10(mean(spike_rates_sessions{a2}(:,LOC1:LOC2),2)./mean(spike_rates_sessions{a2}(:,bs1:bs2),2));
        [p h] = ranksum(tmp1, tmp2)
        effect_size_compare(a1,a2) = p;
    end
end
        

figure;
subplot(2,1,1);
bar(1:4, [FR_awake(:,1) FR_LOC(:,1)]);
hold on;
errorbar([1:4]-0.15, FR_awake(:,1), FR_awake(:,2), 'linestyle', 'none');
errorbar([1:4]+0.15, FR_LOC(:,1), FR_LOC(:,2), 'linestyle', 'none');
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'});
ylabel('Firing Rate (spikes/sec)')
xlim([0.5 4.5])
subplot(2,1,2);
bar(1:4, effect_size);
hold on; errorbar(1:4, effect_size, effect_size_sem, 'linestyle', 'none');
set(gca, 'xticklabel', {'PFC', '8A', 'PPC', 'STG'});
ylabel('Effect Size (dB change from baseline)')
xlim([0.5 4.5])







