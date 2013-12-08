function SegmentPSDSyncAve(PCG,PCGS1BeginIndex, PCGS1EndIndex, Fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% --------------------------------------- Step 1: Find PSDs of each segment
SegmentPSD = zeros(length(PCGS1BeginIndex), 2048);

figure;
for SegmentCounter = 1: length(PCGS1BeginIndex);

    CurrentPCGSegment = PCG( PCGS1BeginIndex(SegmentCounter) : PCGS1EndIndex(SegmentCounter) );
    
    % Subtract the mean of each signal segment
    CurrentPCGSegment = CurrentPCGSegment - mean(CurrentPCGSegment);            
    
    % Compute the PSD of each segment, given as the square of the absolute value of the FT 
    % of the signal segment, divided by the number of samples in the FT array. 
    %figure;
    SegmentPSD(SegmentCounter,:) = FindPSD(CurrentPCGSegment, 2048, 1000);       
end

% --------------------------------------- Step 2: Sync Ave. of segment PSDs
%Obtain the averaged PSD of the systolic portions of each PCG signal
%using as many cardiac cycles as possible in a synchronized averaging procedure
SyncAveragePSD = sum(SegmentPSD,1) / length(PCGS1BeginIndex);

SyncAverageAxis = (0 : (ceil(length(SyncAveragePSD)) - 1 )) * ((Fs) / length(SyncAveragePSD));
SyncAveragePSD = SyncAveragePSD(1:length(SyncAverageAxis));

figure;
plot(SyncAverageAxis, SyncAveragePSD);
title('Sync Average PSD');
xlabel('Freq (Hz)'); ylabel('Magnitude'); axis auto;

% --------------------------------------------- Step 3: Calc Mean Frequency
% Compute the mean frequency of the averaged PSD for each case of the textbook). 
% Use only one half of the PSD (for positive frequency only) to compute the mean frequency.

N = size(SyncAveragePSD);
N = N(2)

RunningSum = 0;
for SyncAverageCounter = 1 : N
    
    RunningSum = RunningSum + abs(SyncAveragePSD(SyncAverageCounter));
    
end
TotalEnergy = RunningSum / N

K = 0:round(N / 2) - 1;
RunningSum = 0;
for SyncAverageCounter = 1 : round(N / 2)

    RunningSum = RunningSum + K(SyncAverageCounter) * SyncAveragePSD(SyncAverageCounter);
    
end
MeanFrequency = ((2 * Fs) / (N * TotalEnergy)) * RunningSum

end

