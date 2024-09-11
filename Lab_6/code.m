% Sixth tutorial.
close all; clear all; clc;

load('EMG.mat'); % Load the EMG signal
L = length(EMG); % Duration of the signal in samples
 
Fs = 2500; % Sample frequency in Hz
t_ax = (0:L-1)/Fs; % Time axis of the signal in seconds
% 
% figure(1), plot(t_ax,EMG);
% title('Intramuscular EMG signal');
% xlabel('Time [s]');
% ylabel('AU');

F_EMG = fftshift(fft(EMG)); % Find the DFT of the EMG 
f_ax = (-L/2:L/2-1)*Fs/L;
figure(2), plot(f_ax,abs(F_EMG));
xlabel('Frequency (Hz)');
ylabel('Magnitude(AU)');
title('Spectrum of the signal before the filtering')

%% Create band-stop filters to remove 60 Hz and harmonics
N = 3; % Filter order
band1 = [58 62];
[B1,A1] = butter(N,band1/(Fs/2),'stop'); % Generate filter coefficients
band2 = [118, 122];
[B2, A2] = butter(N, band2/(Fs/2),'stop'); % added
band3 = [178, 182];
[B3, A3] = butter(N, band3/(Fs/2), 'stop'); %added

% Analyze the properties of the filter
H1 = tf(B1,A1,1/Fs,'variable','z^-1'); % Create transfer function object
[z,p,k] = tf2zp(B1,A1); % Calculate zeros and poles
% figure(3), freqz(B1,A1);
% title('Analysis in frequency domain of the filter that remove 60Hz-band frequency ');
figure(4), subplot(3,1,1),zplane(B1,A1);
title('z-plane to represent zeros and poles of the filter that remove 60Hz-band frequency ');
 
H2 = tf(B2,A2,1/Fs,'variable','z^-1'); 
[z,p,k] = tf2zp(B2,A2); % Calculate zeros and poles
% figure(5), freqz(B2,A2);
% title('Analysis in frequency domain of the filter that remove 120Hz-band frequency ');
figure(4), subplot(3,1,2), zplane(B2,A2);
title('z-plane to represent zeros and poles of the filter that remove 120Hz-band frequency ');

H3 = tf(B3, A3, 1/Fs,'variable','z^-1' );
[z,p,k] = tf2zp(B3,A3); % Calculate zeros and poles
% figure(7), freqz(B3,A3);
% title('Analysis in frequency domain of the filter that remove 180Hz-band frequency ');
figure(4), subplot(3,1,3), zplane(B3,A3);
title('z-plane to represent zeros and poles of the filter that remove 180Hz-band frequency ');

%% Filter the signal
% figure(9), hold on, plot(t_ax,EMG,'b'); %original signal
 
EMG_f = filter(B1,A1,EMG);
% figure(9), hold on, plot(t_ax,EMG_f,'r');

% figure(10); hold on, plot(t_ax,EMG,'b');
EMG_f2 = filter(B2,A2,EMG);
% figure(10), hold on, plot(t_ax,EMG_f2,'r');

% figure(11);hold on, plot(t_ax,EMG,'b');
EMG_f3 = filter(B3, A3, EMG);
% figure(11), hold on, plot(t_ax, EMG_f3,'r');

%% DFT of the filtered signal
figure(12);

F_EMG_f = fftshift(fft(EMG_f));
subplot(3,1,1)
plot(f_ax, abs(F_EMG_f),'r');

F_EMG_f2 = fftshift(fft(EMG_f2));
subplot(3,1,2)
plot(f_ax, abs(F_EMG_f2), 'b');

F_EMG_f3 = fftshift(fft(EMG_f3));
subplot(3,1,3);
plot(f_ax, abs(F_EMG_f3), 'g');

%% cascade
figure(13);
subplot(4,1,1);
plot(f_ax,abs(F_EMG),'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude(AU)');

subplot(4,1,2);
hold on;
plot(f_ax,abs(F_EMG),'r','DisplayName','Original');
plot(f_ax,abs(fftshift(fft(EMG_f))),'b','DisplayName','60Hz Filter Applied');
hold off;
xlabel('Frequency (Hz)');
ylabel('Magnitude(AU)');
legend;

EMG_cascade_1 = filter(B2,A2,EMG_f);
subplot(4,1,3);
hold on;
plot(f_ax,abs(F_EMG),'r','DisplayName','Original');
plot(f_ax,abs(fftshift(fft(EMG_cascade_1))),'b','DisplayName','Cascade 2 Filters');
xlabel('Frequency (Hz)');
ylabel('Magnitude(AU)');
hold off;
legend;

EMG_cascade = filter(B3,A3, EMG_cascade_1);
F_EMG_cascade = fftshift(fft(EMG_cascade));
subplot(4,1,4);
hold on;
plot(f_ax,abs(F_EMG),'r','DisplayName','Original');
plot(f_ax,abs(F_EMG_cascade),'b','DisplayName','Cascade all Filters');
hold off;
xlabel('Frequency (Hz)');
ylabel('Magnitude(AU)');
legend;