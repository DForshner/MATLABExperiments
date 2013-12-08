% Peak detection
% Generates a 10 second sample of signal x(t) at fs=100 Hz.
% which is passed to the function Peak_Detect which returns a vector 
% containing the sample numbers of the peaks & valleys of the signal.

clear all % Clear all Varables
close all % Close all windows

% Generate a sampled version of the signal at the sampling frequency 
% fs = 100 Hz for the duration of 0 − 10 seconds.
Frequency = 100;
Period = 1 / Frequency;
t = 0 : Period : 10; % Ten seconds in increments of Period

x = 3*sin(2*pi*5*t) + 5*cos(2*pi*7*t);  % Create sampled version of X(t)

PEAKS = Peak_Detect(x, 2, 1000);  % Call function to find the peaks & valleys

% Create Plot of sampled signal with peaks & valleys highlighted
figure('Name','Signal with peaks & valleys marked');
plot(t, x, 'k-'); hold on; % Plot the signal x
plot(t(PEAKS), x(PEAKS), 'r*'); hold off; % Highlight Peaks & Valleys
title('Signal with peaks & valleys marked');
xlabel('time (sec)');
ylabel('x(t) = 3*sin(2*pi*5*t) + 5*cos(2*pi*7*t)');