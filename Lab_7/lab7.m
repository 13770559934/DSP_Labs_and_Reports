% Seventh tutorial.
close all; clear all; clc
 
load('EEG.mat') % Load the EEG signal
fs = 200; % Sample frequency in Hz
L = length(EEG); % Duration of the signal in samples
 
winStep = 0.025; % window size increase in seconds
winSizeVector = floor([1:(L)/(winStep*fs)]*winStep*fs); 
% Vector containing window sizes ranging in 25ms step increases
 
% Iterate over different window lengths to compute the estimations of mean
% and variance
for iter_winSize = 1:numel(winSizeVector)
    % calculate the remaining values for the same window size
    for iterSig=1:floor((length(EEG))/winSizeVector(iter_winSize))
        start_idx=(iterSig-1)*winSizeVector(iter_winSize)+1; 
        %starting sample index of the extracted window
        end_idx = iterSig*winSizeVector(iter_winSize); 
        %ending sample index of the extracted window
        V_fixWin(iterSig) = var(EEG(start_idx:end_idx));
        % Calculate variance for the given window length
        M_fixWin(iterSig) = mean(EEG(start_idx:end_idx));
        % Calculate mean for the given window length
    end
    V{iter_winSize,:} = V_fixWin;
    % store variance for the given window length
    M{iter_winSize,:} = M_fixWin;
    % store mean for the given window length
    
    V_avg(iter_winSize) = mean(V_fixWin);% calculate the mean value of the variance for the given window length
    M_avg(iter_winSize) = mean(M_fixWin);% calculate the mean value of the mean for the given window length
        
    V_fixWin=[];
    M_fixWin=[];
end
%%
n = 1;
for i = [0.5,1,1.5]
winSizeIdx=find(winSizeVector == i*fs);
x_axis=[1:length(V{winSizeIdx,:})].*(winSizeVector(winSizeIdx)/fs);
figure(1); 
subplot(3,1,n);
plot(x_axis,V{winSizeIdx,:},'k','LineWidth',2);
xlabel('Time (s)'), ylabel('[AU]');
title(['Variance for the fixed window length = ' num2str(i) 's']);
disp(var(V{winSizeIdx,:}))

figure(2);
subplot(3,1,n);
plot(x_axis, M{winSizeIdx,:},'r','LineWidth',2);
xlabel('Time (s)'), ylabel('[AU]');
title(['Mean for the fixed window length = ' num2str(i) 's']);
n = n+1;
end


x_ax_length = [1:(L)/(winStep*fs)]*winStep; %build the x axis for the plot
% Plot the variance with respect to the window length
figure(3)
subplot(2,1,1)
plot(x_ax_length,V_avg,'k','LineWidth',2),hold on
xlabel('Duraition of the windows (s)'), ylabel('[AU]')
title('Values of variance depending on the window size')


x_ax_length = [1:(L)/(winStep*fs)]*winStep; %build the x axis for the plot
% Plot the variance with respect to the window length
subplot(2,1,2);
plot(x_ax_length,M_avg,'k','LineWidth',2),hold on
xlabel('Duraition of the windows (s)'), ylabel('[AU]')
title('Values of mean depending on the window size')

 
figure
histogram(EEG,'Normalization','probability');
xlabel('EEG(AU)'), ylabel('Estimate of the probability(AU)')
title('PDF of EEG')
