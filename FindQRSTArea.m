function QRSTArea = FindQRSTArea(ECGVoltage, ECGTimeAxis, fs, StartingSample, EndingSample)
% FindQRS - This function calculates QRSTArea
% This function computes QRSTArea given the current peak sample number, the original signal and fs.
%
% INPUTS:
% OUTPUTS:

% Calculate some signal parameters
T = 1/fs; % period

ECGSection = ECGVoltage(StartingSample : EndingSample);

% Subtract the mean from the section of signal
ECGSection = ECGSection - mean(ECGSection);

% Rectify the section of signal
ECGSection = abs(ECGSection);

% Calculate the area
QRSTArea = sum(ECGSection) * T;




