% Eighth tutorial.
close all; clear all; clc
 
load('spike_neural.mat') % Load the neural_sig signal
L = length(neural_sig); % Duration of the signal in samples
fs = 10240; % Sample frequency in Hz
WinSize_1 = [0.2, 0.5, 2]; % Window size in seconds 
WinSize = round(WinSize_1.*fs); % Window size in samples
 
f_ax = (-pi:2*pi/fs:pi-2*pi/fs)./(2*pi).*fs; % Frequency axis in Hz
variance_periodogram_estimate = [];
mean_periodogram = [];
bias_values = [];
figure(1);
for uu = 1 : length(WinSize)
    
    % Estimate of the periodogram
    window = rectwin(WinSize(uu))';
    for n = 1:35 %For each window length, estimate the periodogram for the first 35 signal segments
        wind_signal=neural_sig((n-1)*WinSize(uu)+(1:WinSize(uu))).*window;
        Segm_spect{uu}(n,:) = fftshift(abs(fft(wind_signal,fs)).^2)./WinSize(uu);% calculating PSD
    end    
    
    
    mean_periodogram = [mean_periodogram mean(Segm_spect{uu},1)'];
    figure(1);
    subplot(3,1,uu)
    plot(f_ax,mean_periodogram(:,uu),"DisplayName",sprintf('Window length=%d',WinSize_1(uu)));
    xlim([-50,50])
    title(['Mean Periodogram of Rectangular Window, Length =' num2str(WinSize_1(uu)) ' s']);
    xlabel('Frequency(Hz)')
    ylabel('PSD(AU)')

    freq_range = [-5, 5]; % Hz
    freq_indices = find(f_ax >= freq_range(1) & f_ax <= freq_range(2));
    psd_area = trapz(f_ax(freq_indices),mean_periodogram(freq_indices,uu));
    peak_value = max(mean_periodogram(freq_indices,uu));
    bias_values(uu) = psd_area / peak_value;
end



%%

WinSize = [0.2:0.1:2]; % Window size in seconds 
WinSize_samples = round(WinSize.*fs); % Window size in samples

% 定义窗口类型
windows = {'rectwin', 'hann', 'hamming'};
window_labels_1 = {'Rectangular', 'Hanning', 'Hamming'};
window_labels = arrayfun(@(x) sprintf('%.1fs', x), WinSize, 'UniformOutput', false); % 窗口标签

% 初始化变量以存储方差数据
variances_rect = [];
variances_hann = [];
variances_hamm = [];

% 计算每种窗口类型和大小的方差
for w = 1:length(windows)
    for uu = 1:length(WinSize_samples)
        window = feval(windows{w}, WinSize_samples(uu))';
        Segm_spect{uu} = zeros(35, length(f_ax)); % 假设有35个信号段
        for n = 1:35
            % do FFT?
            wind_signal = neural_sig((n-1)*WinSize_samples(uu)+(1:WinSize_samples(uu))).*window;
            Segm_spect{uu}(n,:) = fftshift(abs(fft(wind_signal, length(f_ax))).^2)./WinSize_samples(uu);
            % calculate PSD
            
        end
        variances{w}(:,uu) = var(Segm_spect{uu})';
    end

end

% 绘制箱形图
figure(2);

% 为每种窗口类型绘制箱形图
for w = 1:length(windows)
    subplot(3, 1, w);
    boxplot(variances{w}, 'Labels', window_labels);
    title(['Variance of Estimation of the Periodogram for ' window_labels_1{w} ' Window']);
    ylabel('Variance Estimation');
    if w == length(windows)
        xlabel('Window Size');
    end
end

% 设置背景为白色
set(gcf, 'Color', 'w');