% Fifth Tutorial
clear all; close all; clc;
 
%% prepare
load 'emg'; % Load the EMG signal
fs = 1600; % Sampling frequency
L = length(emg); % Duration of the signal in samples
time_ax=[0:1/fs:(L-1)/fs]; % Time axis of the signal in seconds
 
fc = 5; % Filter cut-off frequency in (Hz)
wc = 2*fc/fs*pi; % Normalized cut-off frequency
 
%% Generate Sinc function with a given sampling rate
t = -floor(L):floor(L); % Form the time axis of the Sinc function (theoretically it can be infinetly long)
sinc_func=(wc/pi)*sinc((2*fc/fs)*t);
 
figure(1); subplot(2,1,1), plot(t, sinc_func);
title('Sinc function');
xlabel('Discrete samples (n)');
ylabel('AU');
 
f_snc=fftshift(fft(sinc_func)); % Find the DFT of the sinc function 
f_ax =((-pi+pi/length(sinc_func):2*pi/length(sinc_func):pi-pi/length(sinc_func))./pi) .* (fs/2); % Frequency axis for the DFT of sinc function

figure(1);subplot(2,1,2), plot(f_ax,abs(f_snc));
title('Magnitude of discrete time Fourier transform of the Sinc');
xlabel('Frequency (Hz)');
ylabel('Magnitude (AU)');
xlim([-20, 20])
 
%% Create the FIR filter by truncating the very long sinc function and by using
% one of the following windows: hanning (hann) and rectangular (rectwin);
filter_durations = [400,1000];
for i = [1 2]
    FIR_duration=filter_durations(i);
    truncation_section=floor(length(sinc_func)/2-FIR_duration/2):floor(length(sinc_func)/2+FIR_duration/2);
    fir_filt = sinc_func(truncation_section).*rectwin(length(truncation_section))';
    
    
    figure(4+i);freqz(fir_filt,1,8192);
    
    % Create moving avarege filter with the lenght of MA_coef_num
    MA_coef_num=filter_durations(i);
    MA = ones(1,MA_coef_num)/MA_coef_num; % Impulse response of the moving average filter (see slide 30)
    figure(6+i);freqz(MA); 
    
    emg = abs(emg);
    % Filter the rectified EMG using FIR filter
    env_FIR = conv(emg, fir_filt,'same');
    figure(2), subplot(2,1,i), plot(time_ax, env_FIR);
    title(sprintf('EMG signal filtered with FIR filter, length = %d', FIR_duration));
    xlabel('Time (s)');
    ylabel('EMG (AU)');
    
    % Filter the rectified EMG using moving avarage filter
    env_MA = conv(emg,MA,'same');
    figure(3), subplot(2,1,i), plot(time_ax, env_MA);
    title(sprintf('EMG signal filtered with MA filter, length = %d', FIR_duration));
    xlabel('Time (s)');
    ylabel('EMG (AU)');
end
FIR_duration=400;
truncation_section=floor(length(sinc_func)/2-FIR_duration/2):floor(length(sinc_func)/2+FIR_duration/2);
fir_filt = sinc_func(truncation_section).*rectwin(length(truncation_section))';
figure(4);
subplot(2,1,1);
env_FIR = conv(emg, fir_filt,'same');
plot(time_ax, env_FIR);
title(sprintf('EMG signal filtered with FIR filter, length = %d, rectangular window', FIR_duration));
xlabel('Time (s)');
ylabel('EMG (AU)');
fir_filt = sinc_func(truncation_section).*hann(length(truncation_section))';
env_FIR = conv(emg, fir_filt,'same');
subplot(2,1,2);
plot(time_ax, env_FIR);
title(sprintf('EMG signal filtered with FIR filter, length = %d, hanning window', FIR_duration));
xlabel('Time (s)');
ylabel('EMG (AU)');
