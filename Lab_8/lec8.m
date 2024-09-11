% Eighth tutorial.
close all; clear all; clc
 
load('spike_neural.mat') % Load the neural_sig signal
L = length(neural_sig); % Duration of the signal in samples
fs = 10240; % Sample frequency in Hz
WinSize = [0.2:0.1:2]; % Window size in seconds 
WinSize = round(WinSize.*fs); % Window size in samples
 
f_ax = (-pi:2*pi/fs:pi-2*pi/fs)./(2*pi).*fs; % Frequency axis in Hz
variance_periodogram_estimate = [];
for uu = 1 : length(WinSize)
    
    % Estimate of the periodogram
    window = rectwin(WinSize(uu))';
    for n = 1:35 %For each window length, estimate the periodogram for the first 35 signal segments
        wind_signal=neural_sig((n-1)*WinSize(uu)+(1:WinSize(uu))).*window;
        Segm_spect{uu}(n,:) = fftshift(abs(fft(wind_signal,fs)).^2)./WinSize(uu);
    end    
    
    %Variance of the periodogram estimate for each window size
    variance_periodogram_estimate = [variance_periodogram_estimate var(Segm_spect{uu})'];
 
end

% Define specific window sizes in seconds
specific_windows_sec = [0.2, 0.5, 2]; 
specific_windows_samples = round(specific_windows_sec .* fs);

% Frequency range for plotting
freq_range = -50:50;

% Calculate the indices for the desired frequency range
freq_indices = find(f_ax >= -50 & f_ax <= 50);

% Initialize figure
figure;

% Plotting mean periodogram for each specific window size
for w = 1:length(specific_windows_sec)
    window_size = specific_windows_samples(w);
    % Find index of the window size in WinSize array
    idx = find(WinSize == window_size);
    % Calculate mean periodogram
    mean_per = mean(Segm_spect{idx}, 1);
    % Select the desired frequency range from mean periodogram
    mean_per_freq_range = mean_per(freq_indices);
    % Scale PSD for plotting
    mean_per_freq_range_scaled = mean_per_freq_range / 1e3;
    % Plot
    subplot(length(specific_windows_sec), 1, w);
    plot(f_ax(freq_indices), mean_per_freq_range_scaled);
    xlim([-50 50]);
    ylim([min(mean_per_freq_range_scaled), max(mean_per_freq_range_scaled)]); % Adjust y-limits based on your data
    xlabel('Frequency / Hz');
    ylabel('PSD');
    title(['Mean Periodogram for Rectangular Window Size at ' num2str(specific_windows_sec(w)) ' s']);
end


% Specific window sizes in seconds and their corresponding samples
specific_windows_sec = [0.2, 0.5, 2];
specific_windows_samples = round(specific_windows_sec .* fs);

% Frequency range for bias computation
freq_range_bias = -5:5;

% Pre-allocate table for bias reporting
bias_table = array2table(zeros(length(specific_windows_sec), 1), ...
                         'VariableNames', {'Bias'}, ...
                         'RowNames', {'0.2s', '0.5s', '2s'});

% Loop through each window size to compute bias
for w = 1:length(specific_windows_sec)
    window_size_samples = specific_windows_samples(w);
    
    % Find index of the window size in WinSize array
    idx = find(WinSize == window_size_samples);
    
    % Calculate mean periodogram for the window size
    mean_per = mean(Segm_spect{idx}, 1);
    
    % Identify peak value within the specified frequency range
    freq_indices_bias = find(f_ax >= -5 & f_ax <= 5);
    peak_value = max(mean_per(freq_indices_bias));
    
    % Compute the area under the curve around the peak
    area_around_peak = trapz(f_ax(freq_indices_bias), mean_per(freq_indices_bias));
    
    % Normalize this area by the peak value
    normalized_area = area_around_peak / peak_value;
    
    % Store the normalized area in the table
    bias_table.Bias(w) = normalized_area;
end

% Display the table
disp(bias_table)
