function QRSBeginIndex = FindQRSBegin(DownSampledECG, DownSampledTimeAxis, QRSPeakIndex, Fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[b,a] = butter(2, (1.5/ (Fs/2)), 'high'); % Cut out anything below 1.5 Hz to optimize search forward routine
NoBiasDownSampledECG = filter(b, a, DownSampledECG);

QRSBeginIndex = []; % Index of QRS Start samples

for QRSPeakIndexCounter=1:length(QRSPeakIndex)

    CurrentSample = QRSPeakIndex(QRSPeakIndexCounter);
    
    while (NoBiasDownSampledECG(CurrentSample) > (0.20 * NoBiasDownSampledECG(QRSPeakIndex(QRSPeakIndexCounter))) && (CurrentSample > 0))
    
        CurrentSample = CurrentSample - 1;
    
    end
   
    QRSBeginIndex(QRSPeakIndexCounter) = CurrentSample;
    
end

QRSStartPlot = figure('Name','Start of QRS - ECG signal'); % Create a new figure
plot(DownSampledTimeAxis, DownSampledECG); 
hold on; plot(DownSampledTimeAxis(QRSBeginIndex), DownSampledECG(QRSBeginIndex), 'r*'); hold off;
title('QRS start locations - ECG signal (w/ starting points indicated by red "*" ');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis auto;

end

