% Third tutorial
clear all
close all
clc
 
% Signal loading
load('EEG.mat');
 
% Sampling frequency
fsamp = 512;

% Select duration to analyze
Durations = [1,5,15]; % Duration in seconds (max 15 seconds) in sec??
figure(1); 
figure(2);
sgtitle('Power Spectral Density of EEG')

for Duration = Durations
    load('EEG.mat');
    Duration = round(Duration*fsamp);
    EEG = EEG(1:Duration); 
     
    % Signal duration in samples
    L = length(EEG);
     
    % Plot the EEG signal
    time_ax = [0:L-1]./fsamp;
    figure(1)
    plot(time_ax, EEG);
    xlabel('Time (s)')
    title(['EEG signal'])
    ylabel('Amplitude EEG (Arbitrary Units)')
     
    % Compute DFT 
    X1 = fft( EEG - mean(EEG) );
     
    % Compute PSD (power spectral density)and adjust
    % frequency axis according to Matlab notation
    PSD1 = fftshift(abs(X1).^2)/L;
     
    % Build the frequency axis in radiants
    freq_a_rad = [-pi+pi/L:2*pi/L:pi-pi/L];
    % Convert the frequency axis in Hz
    freq_a_Hz = freq_a_rad./(2*pi).*fsamp;
    
    figure(2), subplot(2,3,find(Duration/fsamp == Durations)), plot(freq_a_Hz,PSD1);
    xlabel('Frequency (Hz)');
    ylabel('PSD (Arbitrary Units)');
    title(sprintf('Signal interval = %d(s)', Duration/fsamp));
    subplot(2,3,find(Duration/fsamp == Durations) + 3)
    plot(freq_a_rad,PSD1);
    xlabel('Frequency (radians)')
    ylabel('PSD (Arbitrary Units)')

    %% Compute the percentage of power in different subbands
    %freq_a_Hz = round(freq_a_Hz .* 10) ./ 10;
    halfDuration = Duration/2;
    total_power = sum(PSD1(halfDuration+1 : end));
    indexs = [[0.5,4];[4,8];[8, 13];[13,30];[30,42]];
    percentages = [];

    for i = 1:size(indexs,1)
        index = find(freq_a_Hz >= indexs(i,1) & freq_a_Hz <= indexs(i,2));
        PSDs = PSD1(index);
        segement_power = sum(PSDs);
        percentage = (segement_power/total_power) * 100;
        percentages = [percentages, percentage];
    end
    disp(Duration/fsamp)
    disp(percentages);
end



%% Compute the percentage of power in different subbands
freq_a_Hz = round(freq_a_Hz .* 10) ./ 10;
halfDuration = Duration/2;
total_power = sum(PSD1(halfDuration+1 : end));
indexs = [[0.5,4];[4,8];[8, 13];[13,30];[30,42]];
percentages = [];

for i = 1:size(indexs,1)
    index = find(freq_a_Hz >= indexs(i,1) & freq_a_Hz <= indexs(i,2));
    PSDs = PSD1(index);
    segement_power = sum(PSDs);
    percentage = (segement_power/total_power) * 100;
    percentages = [percentages, percentage];
end