figure(1);
sgtitle('MSE against Delay for different M factor');
for M = [1 2 4 8]
    n_step = 100; % Number of steps in the loop for optimal alignment
    stepSize = 0.2; % Size of the step of non-integer delay applied to one of the signals for alignment
    ied=24; % Interelectrode distance of the recordings in mm
    Fs = 2048; % Recording sampling frequency
    Ts = 1/Fs; % Sampling interval
    load('Signals_EMG.mat'); % Loading the recorded EMGs (two channels)

    channel1 = channel1(1:M:end); % Downsampling first channel
    channel2 = channel2(1:M:end); % Downsampling second channel
    
    Fs=Fs/M;
    Ts=Ts*M;
    % End downsampling

    timeAxis=[1:length(channel1)].*Ts.*1000; % Definition of time axis in ms
    freqAxis=fftshift([-0.5:1/(length(channel1)):0.5-1/(length(channel1))]); % Definition of discrete frequency axis
 
    channel1_ft = fft(channel1); % Fourier transform of the first channel
    
    for uu = 1 : n_step
        channel1_dt = (channel1_ft).*exp(-i*2*pi*stepSize*uu*freqAxis); % complex exponential multiplication (delay in frequency domain)
        channel1_dt = real(ifft((channel1_dt))); % inverse transform to go back to the time domain 
        MSE_vect(uu)= sum((channel1_dt - channel2).^ 2)./sum(channel2.^ 2).*100; % normalized mean square error between aligned signals
        delay(uu) = (stepSize*uu); % Imposed delay in samples
    end

    [MSEopt, optDelay] = min(MSE_vect);
 
    % Delay and conduction velocity estimate
    fprintf('The optimal delay is %2.2f ms \n',delay(optDelay)*Ts*1000);
    fprintf('The estimated conduction velocity is %2.2f m/s \n',ied/(delay(optDelay)*Ts*1000));
    fprintf('Optimal MSE between signals: %2.2f %%\n',MSEopt);
    
    subplot(2,2,((log(M) / log(2)+1)));
    plot( delay .* (Ts * 1000), MSE_vect);

    title(sprintf("M = %d", M))
    xlabel('Delay(ms)');
    ylabel('Estimation Error(%)');
    clear all;
end