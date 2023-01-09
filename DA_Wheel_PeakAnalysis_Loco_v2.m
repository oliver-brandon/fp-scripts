%DA_Wheel_PeakAnalysis_Loco.m created by Brandon L. Oliver, M.A. for use by the
%lab of Dr. Natalie Zlebnik

%Input = Locomotor DA recording tanks
%Output = table called "peak_analysis" that includes number of peaks,
%peaks/s, max amp, and avg amp for DLS and NAc signals

%Instructions: Point MATLAB at the desired TDT tank folder (data), set
%session time in seconds (session_duration), and click run. Data is output
%as a table called "peak_analysis."

%Find peaks in the 465 signal that are above the median absolute deviation
%(MAD)for the signal. Users can also specifiy a minimum peak distance threshold. 
%MAD1 = MAD value for the 465A and MAD2 = MAD value for the 465C

%NEW: This script now supports batch analysis of tanks in a folder. Simply
%run the script and choose the folder containing the tanks you wish to
%analyze.

%House Keeping%
clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
warning

%Choose batch or single file analysis%
tanks2analyze = 1; %1 = batch analyze, 2 = single file analysis%
%Choose minimum peak distance%
peakDist = 0.2;
%Choose session duration%
session_duration = 3600; %duration of recording in seconds

if tanks2analyze == 1
    myDir = uigetdir; %gets directory%
    myFiles = dir(myDir); %gets all tanks in directory%
    myFiles = myFiles(~ismember({myFiles.name},{'.','..'}));
    for i = 1:length(myFiles)
        BLOCKPATH = fullfile(myDir, myFiles(i).name);
        data = TDTbin2mat(BLOCKPATH, 'T2', session_duration, 'TYPE', {'streams'});
        
        %Stream Stores%
        DLS_ISOS = 'x405A'; % name of the 405A store
        DLS_DA = 'x465A'; % name of the 465A store
        NAc_ISOS = 'x405C'; % name of the 405C store
        NAc_DA = 'x465C'; % name of the 465C store
        N = 100; %Downsample N times
        %time array used for all streams%
        time = (1:length(data.streams.(DLS_DA).data))/data.streams.(DLS_DA).fs;
        %removes the first (t) seconds where the data is wild due to turning on LEDs%
        t = 20; % time threshold below which we will discard
        ind = find(time>t,1);% find first index of when time crosses threshold
        time = time(ind:end); % reformat vector to only include allowed time
        data.streams.(DLS_DA).data = data.streams.(DLS_DA).data(ind:end);
        data.streams.(DLS_ISOS).data = data.streams.(DLS_ISOS).data(ind:end);
        data.streams.(NAc_DA).data = data.streams.(NAc_DA).data(ind:end);
        data.streams.(NAc_ISOS).data = data.streams.(NAc_ISOS).data(ind:end);
        
        %downsample streams and time array by N times%
        data.streams.(DLS_ISOS).data = downsample(data.streams.(DLS_ISOS).data, N);
        data.streams.(DLS_DA).data = downsample(data.streams.(DLS_DA).data, N);
        data.streams.(NAc_ISOS).data = downsample(data.streams.(NAc_ISOS).data, N);
        data.streams.(NAc_DA).data = downsample(data.streams.(NAc_DA).data, N);
        time = downsample(time, N);
        
        %detrend & dFF%
        %465A%
        bls = polyfit(data.streams.(DLS_ISOS).data,data.streams.(DLS_DA).data,1);
        Y_fit_all = bls(1) .* data.streams.(DLS_ISOS).data + bls(2);
        Y_dF_all = data.streams.(DLS_DA).data - Y_fit_all; %dF (units mV) is not dFF
        dFF = 100*(Y_dF_all)./Y_fit_all;
        std_dFF = std(double(dFF));
        detrend_465A = detrend(dFF);
        detrend_465A = zscore(detrend_465A);
        %465C%
        bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
        Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
        Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF
        dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
        std_dFF2 = std(double(dFF2));
        detrend_465C = detrend(dFF2);
        detrend_465C = zscore(detrend_465C);
        %calculates and plots median absolute deviation for both 465 signals%
        MAD1 = mad(detrend_465A, 1);
        MAD2 = mad(detrend_465C, 1);
        [pks,locs,w,p] = findpeaks(detrend_465A, time, 'MinPeakHeight', MAD1,...
            "MinPeakDistance", peakDist);
        [pks2,locs2,w2,p2] = findpeaks(detrend_465C, time, 'MinPeakHeight', MAD2,...
            "MinPeakDistance", peakDist);
        
        DLS_pks(i,:) = length(pks);
        DLS_pk_min(i,:) = (DLS_pks(i)/session_duration)*60;
        DLS_amp_max(i,:) = max(pks);
        DLS_amp_avg(i,:) = mean(pks);
        NAc_pks(i,:) = length(pks2);
        NAc_pk_min(i,:) = (NAc_pks(i)/session_duration)*60;
        NAc_amp_max(i,:) = max(pks2);
        NAc_amp_avg(i,:) = mean(pks2);

    end
elseif tanks2analyze == 2
    myDir = uigetdir;
    BLOCKPATH = myDir;
    data = TDTbin2mat(BLOCKPATH, 'T2', session_duration, 'TYPE', {'streams'});
    
    %Stream Stores%
    DLS_ISOS = 'x405A'; % name of the 405A store
    DLS_DA = 'x465A'; % name of the 465A store
    NAc_ISOS = 'x405C'; % name of the 405C store
    NAc_DA = 'x465C'; % name of the 465C store
    N = 100; %Downsample N times
    %time array used for all streams%
    time = (1:length(data.streams.(DLS_DA).data))/data.streams.(DLS_DA).fs;
    %removes the first (t) seconds where the data is wild due to turning on LEDs%
    t = 20; % time threshold below which we will discard
    ind = find(time>t,1);% find first index of when time crosses threshold
    time = time(ind:end); % reformat vector to only include allowed time
    data.streams.(DLS_DA).data = data.streams.(DLS_DA).data(ind:end);
    data.streams.(DLS_ISOS).data = data.streams.(DLS_ISOS).data(ind:end);
    data.streams.(NAc_DA).data = data.streams.(NAc_DA).data(ind:end);
    data.streams.(NAc_ISOS).data = data.streams.(NAc_ISOS).data(ind:end);
    
    %downsample streams and time array by N times%
    data.streams.(DLS_ISOS).data = downsample(data.streams.(DLS_ISOS).data, N);
    data.streams.(DLS_DA).data = downsample(data.streams.(DLS_DA).data, N);
    data.streams.(NAc_ISOS).data = downsample(data.streams.(NAc_ISOS).data, N);
    data.streams.(NAc_DA).data = downsample(data.streams.(NAc_DA).data, N);
    time = downsample(time, N);
    
    %detrend & dFF%
    %465A%
    bls = polyfit(data.streams.(DLS_ISOS).data,data.streams.(DLS_DA).data,1);
    Y_fit_all = bls(1) .* data.streams.(DLS_ISOS).data + bls(2);
    Y_dF_all = data.streams.(DLS_DA).data - Y_fit_all; %dF (units mV) is not dFF
    dFF = 100*(Y_dF_all)./Y_fit_all;
    std_dFF = std(double(dFF));
    detrend_465A = detrend(dFF);
    detrend_465A = zscore(detrend_465A);
    %465C%
    bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
    Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
    Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF
    dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
    std_dFF2 = std(double(dFF2));
    detrend_465C = detrend(dFF2);
    detrend_465C = zscore(detrend_465C);
    %calculates and plots median absolute deviation for both 465 signals%
    MAD1 = mad(detrend_465A, 1);
    MAD2 = mad(detrend_465C, 1);
    [pks,locs,w,p] = findpeaks(detrend_465A, time, 'MinPeakHeight', MAD1,...
        "MinPeakDistance", peakDist);
    [pks2,locs2,w2,p2] = findpeaks(detrend_465C, time, 'MinPeakHeight', MAD2,...
        "MinPeakDistance", peakDist);
    
    DLS_pks = length(pks);
    DLS_pk_min = (DLS_pks/session_duration)*60;
    DLS_amp_max = max(pks);
    DLS_amp_avg = mean(pks);
    NAc_pks = length(pks2);
    NAc_pk_min = (NAc_pks/session_duration)*60;
    NAc_amp_max = max(pks2);
    NAc_amp_avg = mean(pks2);
end


loco_peak_analysis = table(DLS_pks, DLS_pk_min, DLS_amp_max, DLS_amp_avg, ...
NAc_pks, NAc_pk_min, NAc_amp_max, NAc_amp_avg, ...
'VariableNames', {'DLS Peaks','DLS Peaks/m','DLS Max Amp', ...
'DLS Avg Amp','NAc Peaks','NAc Peaks/m','Nac Max Amp', ...
'NAc Avg Amp'});

%UITable (figure) that displays "peak_analysis" table%
figure;
uitable('Data',loco_peak_analysis{:,:},'ColumnName',loco_peak_analysis.Properties.VariableNames,...
    'RowName',loco_peak_analysis.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
disp("DONE")