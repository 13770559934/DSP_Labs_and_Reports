%lab_4
clear all; close all; clc;
 
load('ECG.mat') % Load the signal
ECG=ECG-mean(ECG); % Remove mean
fs = 1000; % Sample frequency in Hz
t_ax = (0:length(ECG)-1)/fs; % Time axis of the signal
 
%% Plot the signal
figure(1), subplot(4,1,1),plot(t_ax, ECG );
title('ECG original signal');
xlabel('Time (s)'),ylabel('Amplitude (AU)');
xlim([0 5]);
 
ECG_duration=size(ECG,2); % Duration of the ECG signal in samples
f_ax=[-pi+pi/ECG_duration:2*pi/ECG_duration:pi-pi/ECG_duration]; % Frequency axis for DFT
F_ECG = fftshift(fft(ECG)); % Calculate DFT 
 
%% Plot the DFT of the signal
figure(2), subplot(4,1,1),plot(f_ax,abs(F_ECG));
title('DFT of the ECG original signal');
xlabel('Frequency (rad)');
ylabel('Magnitude (AU)');

MA_coef_nums = [2 ,10,50];

for i = 1:length(MA_coef_nums)
 
    %% Create moving average filter
    MA_coef_num = MA_coef_nums(i);
    MA = ones(1,MA_coef_num)/MA_coef_num; % Impulse response of the moving average filter (see slide 30)
    ECG_filt = conv(ECG,MA,'same');
    
    % Plot the filtered ECG signal 
    figure(1); subplot(4,1,i+1), plot(t_ax, ECG_filt, 'r');
    title(sprintf('Filtered ECG signal with MA filter(length = %d)', MA_coef_num));
    xlabel('Time (s)'),ylabel('Amplitude (AU)');
    xlim([0 5]);
     
    F_ECG_filt = fftshift(fft(ECG_filt));% Calculate DFT of the filtered ECG
     
    % Plot the Fourier transform of the filtered signal
    figure(2), subplot(4,1,i+1),plot(f_ax,abs(F_ECG_filt),'r');
    title(sprintf('DFT of the filtered ECG signal(filter length = %d)',MA_coef_num));
    xlabel('Frequency (rad)');
    ylabel('Magnitude (AU)');
     
    % Create system function of the z-transform of the MA filter
    H = tf(MA,1,1/fs,'variable','z^-1');
     
    % Evaluate the magnitude and the phase of the filter with respect to the normalized frequency 
    figure(5),subplot(3,1,i) ,freqz(MA,1);
    % Go from coefficients to zeroes and poles of the moving average filters
    [MA_zeros,MA_poles] = tf2zpk(H.Numerator{1,1},(H.Denominator{1,1}));
    % Z-plane of the moving average filter
    figure(6),subplot(3,1,i); zplane(MA_zeros,MA_poles);
    title(sprintf('Filter coefficients represented in the z-plane (length=%d)',MA_coef_num));

end

figure(6), sgtitle('Filter coefficients represented in the z-plane with different length');