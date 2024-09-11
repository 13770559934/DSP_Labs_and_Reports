close all;
clear all;
fclose('all');
 
load('Signals_EMG.mat');  % Loading the recorded EMGs (two channels)
 
n_step = 100;  % Number of steps in the loop for optimal alignment
stepSize = 0.2;  % Size of the step of non-integer delay applied to one of the signals for alignment
ied = 24;  % Interelectrode distance of the recordings in mm
Fs = 2048;  % Recording sampling frequency
Ts = 1/Fs;  % Sampling interval
 
% Downsampling by an integer
M = 8;  % Downsampling factor
channel1 = channel1(1:M:end);  % Downsampling first channel
channel2 = channel2(1:M:end);  % Downsampling second channel
Fs = Fs/M;
Ts = Ts*M;
 
% End downsampling
 
timeAxis = [1:length(channel1)].*Ts.*1000;  % Definition of time axis in ms
freqAxis = fftshift([-0.5:1/(length(channel1)):0.5-1/(length(channel1))]);  % Definition of discrete frequency axis
 
figure(1);
plot(timeAxis,channel1,'k');
hold on;
plot(timeAxis,channel2-1000,'k');  % Plot of the two recordings after down-sampling
xlabel('Time (ms)');
ylabel('Signals (AU)');
 
channel1_ft = fft(channel1);  % Fourier transform of the first channel
 
figure(2);
for uu = 1 : n_step
    channel1_dt = (channel1_ft).*exp(-i*2*pi*stepSize*uu*freqAxis);  % complex exponential multiplication (delay in frequency domain)
    channel1_dt = real(ifft((channel1_dt)));  % inverse transform to go back to the time domain
    plot(timeAxis,channel1_dt,'r');  % Plot of the time-shifted signal
    hold on;
    plot(timeAxis,channel2,'k');
    MSE_vect(uu) = sum((channel1_dt - channel2).^ 2)./sum(channel2.^ 2).*100;  % normalized mean square error between aligned signals
    delay(uu) = stepSize*uu;  % Imposed delay in samples
end;
xlabel('Time (ms)');
ylabel('Signal amplitude (AU)');
 
% Identification of the optimal delay (minimum mean sqaure error)
[MSEopt, optDelay] = min(MSE_vect);
 
% Plot of optimal alignment
% Applying the optimal delay to channel1
channel1_opt_dt = channel1_ft .* exp(-i*2*pi*delay(optDelay)*freqAxis);
channel1_opt_dt = real(ifft(channel1_opt_dt));
 
figure(3);
plot(timeAxis, channel1_opt_dt);  % Plot of the optimally time-shifted signal
hold on;
plot(timeAxis, channel2);
xlabel('Time (ms)');
ylabel('Signal amplitude (AU)');
title('Optimal Alignment of Signals');
legend('Channel 1 (Optimally Aligned)', 'Channel 2');
 
% Delay and conduction velocity estimate
fprintf('The optimal delay is %2.2f ms \n',delay(optDelay)*Ts*1000);
fprintf('The estimated conduction velocity is %2.2f m/s \n',ied/(delay(optDelay)*Ts*1000));
fprintf('Optimal MSE between signals: %2.2f %%\n',MSEopt);


for M = [1 2 4 8]
    load("Signals_EMG.mat");
    channel1 = channel1(1:M:end);  % Downsampling first channel
    channel2 = channel2(1:M:end);  % Downsampling second channel
    Fs = Fs/M;
    Ts = Ts*M;
    timeAxis = [1:length(channel1)].*Ts.*1000;  % Definition of time axis in ms
    freqAxis = fftshift([-0.5:1/(length(channel1)):0.5-1/(length(channel1))]);  % Definition of discrete frequency axis
    channel1_ft = fft(channel1);  % Compute the Fourier transform for the downsampled channel1
 
    for uu = 1 : n_step
        channel1_dt = (channel1_ft).*exp(-i*2*pi*stepSize*uu*freqAxis);  % complex exponential multiplication (delay in frequency domain)
        channel1_dt = real(ifft((channel1_dt)));  % inverse transform to go back to the time domain
        MSE_vect(uu) = sum((channel1_dt - channel2).^ 2)./sum(channel2.^ 2).*100;  % normalized mean square error between aligned signals
        delay(uu) = stepSize*uu;  % Imposed delay in samples
    end;
    figure(4);
    plot(delay, MSE_vect);
    hold on;
 
    [MSEopt, optDelay] = min(MSE_vect);
 
    % Delay and conduction velocity estimate
    fprintf('The optimal delay is %2.2f ms \n',delay(optDelay)*Ts*1000);
    fprintf('The estimated conduction velocity is %2.2f m/s \n',ied/(delay(optDelay)*Ts*1000));
    fprintf('Optimal MSE between signals: %2.2f %%\n',MSEopt);
 
end;