close all; clear all; clc
 
load('spike_neural(3).mat') % Load the neural_sig signal
fs = 10240; % Sample frequency in Hz
WinSize_1 = [0.1:0.1:1]; % Window size in seconds 
WinSize = round(WinSize_1.*fs); % Window size in samples
OverlapValues = [0 0.50];

Bias_Estimate_Rectangular = zeros(numel(WinSize),numel(OverlapValues));
Bias_Estimate_Hanning = zeros(numel(WinSize),numel(OverlapValues));
Var_Estimate_Rectangular = zeros(numel(WinSize),numel(OverlapValues));
Var_Estimate_Hanning = zeros(numel(WinSize),numel(OverlapValues));

for uu = 1 : length(WinSize)
    % every length
    overlap_count=1;
    for Overlap=OverlapValues % size of overlap in percents of window size
        %Rectangular window type
        [var_estimate,bias] = welch_periodogram(neural_sig, WinSize(uu), Overlap, 'rect', fs );
        Bias_Estimate_Rectangular(uu,overlap_count) = bias;
        Var_Estimate_Rectangular(uu,overlap_count) = median(var_estimate);  
        %Hanning window type
        [var_estimate,bias] = welch_periodogram(neural_sig, WinSize(uu), Overlap, 'hann', fs );
        Bias_Estimate_Hanning(uu,overlap_count) = bias;
        Var_Estimate_Hanning(uu,overlap_count) = median(var_estimate);
        
        overlap_count = overlap_count + 1;
    end
 
end

%Plots for the bias of the estimate for different window types and overlaps over different window sizes
figure(1);
hold on;
plot(WinSize_1,Bias_Estimate_Hanning(:,1),':','DisplayName','Hanning window,overlap=0', 'LineWidth', 2)
plot(WinSize_1,Bias_Estimate_Hanning(:,2),':','DisplayName','Hanning window,overlap=0.5', 'LineWidth', 2)
plot(WinSize_1,Bias_Estimate_Rectangular(:,1),'DisplayName','Rectangular window,overlap=0', 'LineWidth', 2)
plot(WinSize_1,Bias_Estimate_Rectangular(:,2),'DisplayName','Rectangular window,overlap=0.5', 'LineWidth', 2)
xlabel('Window Sizes (s)')
ylabel('Bias (AU)')
title('Bias of Different Window Types and Overlap')
legend
hold off 

% Plots for the variability of the estimate for different window types and overlaps over different window sizes
figure(2); hold on;
plot(WinSize_1,Var_Estimate_Hanning(:,1),':','DisplayName','Hanning window,overlap=0', 'LineWidth', 2)
plot(WinSize_1,Var_Estimate_Hanning(:,2),':','DisplayName','Hanning window,overlap=0.5', 'LineWidth', 2)
plot(WinSize_1,Var_Estimate_Rectangular(:,1),'DisplayName','Rectangular window,overlap=0', 'LineWidth', 2)
plot(WinSize_1,Var_Estimate_Rectangular(:,2),'DisplayName','Rectangular window,overlap=0.5', 'LineWidth', 2)
xlabel('Window Sizes (s)')
ylabel('Variability (AU)')
title(' Varaibility of Different Window Types and Overlap')
legend
hold off 


[var_estimate,bias] = welch_periodogram(neural_sig, 0.1*fs,0,'rect',fs);

function [var_estimate,bias] = welch_periodogram(Data, WindowSize, Overlap, WindowType, fs )
   % Estimate of the periodogram
    L = length(Data); % Duration of the signal in samples
    f_ax = (-pi:2*pi/fs:pi-2*pi/fs)./(2*pi).*fs; % Frequency axis in Hz
 
    if strcmp(WindowType,'rect')
        window = rectwin(WindowSize);
    elseif strcmp(WindowType,'hann')
        window = hann(WindowSize);
    end
    
    n=1;
    while max(round((n-1)*WindowSize*(1-Overlap))+(1:WindowSize))<=L
        wind_signal= Data(round((n-1)*WindowSize*(1-Overlap))+(1:WindowSize)).*window';
        Segm_spect(n,:) = fftshift(abs(fft(wind_signal, length(f_ax))).^2)./WindowSize;
        n=n+1;
    end    
    
    %periodogram = [PLEASE COMPLETE];
    var_estimate = median(var(Segm_spect));

    freq_range = [-5, 5]; % Hz
    freq_indices = find(f_ax >= freq_range(1) & f_ax <= freq_range(2));
    psd_area = trapz(f_ax(freq_indices),mean(Segm_spect(:,freq_indices)));
    peak_value = max(mean(Segm_spect));
    bias = psd_area / peak_value;
    
    %figure; plot(f_ax,periodogram); 
end
