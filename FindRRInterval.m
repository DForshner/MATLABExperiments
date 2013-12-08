function RRInterval = FindRRInterval(PeakIndex, CurrentPeak, fs)
% FindQRS - This function takes an ECG signal and returns QRS peak locations
% This function is broken up into 9 seperate stages.
%
% INPUTS: ECGVoltage (ECG voltage Signal), ECGTimeAxis (ECG Time axis), fs(Sampling frequency)
% OUTPUTS: PeakIndex (An index containing the locations of QRS peaks)

% Calculate some signal parameters
T = 1/fs; % period

% Average RR interval and standard deviation of RR intervals of each signal (in ms).

if CurrentPeak == 1
    SamplesBetweenPeaks = 0;
else
    SamplesBetweenPeaks = PeakIndex(CurrentPeak) - PeakIndex(CurrentPeak - 1);   
end

RRInterval = SamplesBetweenPeaks * T; % Calculate the average RR Interval in sec