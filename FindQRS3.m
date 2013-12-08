function QRSIndex = FindQRS(ECGVoltage, ECGTimeAxis, fs)
% FindQRS - This function takes an ECG signal and returns QRS peak locations
% This function is broken up into 9 seperate stages.
%
% INPUTS: ECGVoltage (ECG voltage Signal), ECGTimeAxis (ECG Time axis), fs(Sampling frequency)
% OUTPUTS: PeakIndex (An index containing the locations of QRS peaks)

% Calculate some signal parameters
N = length(ECGVoltage); % number of samples
T = 1/fs; % period

% ----------------------------------------------------------------- Stage 1
% This stage impliments a Low Pass Filter
% Input: ECGVoltage
% Output: ECGLPFilteredVoltage

% Derive low pass (LP) Filter coefficients
% H(z) = (1/32) [1 + z^-6]^2 / [1 + z^-1]^2
% Y(z)[1+z^-1]^2 = (1/32) X(z) [[1 + z^-6]^2 
% Y(z)[1 - 2z^-1 + z^-2] = (1/32) X(z) [1 - 2z^-6 + z^-12]
% y(n) - 2y(n-1) + y(n-2) = (1/32)x(n) - (1/32)2x(n-6) + (1/32)x(n-12)
% put in form a(1)y(n) + a(2)y(n-1) = b(1)x(n) + b(2)x(n-1) + b(7)x(n-6) + b(13)x(n-12)
LPa = [1 -2 1];
LPb = zeros(1,13);
LPb(1) = (1/32); LPb(7) = -(2/32); LPb(13) = (1/32);

% Filter the ECG signal
ECGLPFilteredVoltage = filter(LPb, LPa, ECGVoltage);

% ----------------------------------------------------------------- Stage 2
% This stage impliments a High Pass Filter
% Input: ECGLPFilteredVoltage 
% Output: ECGHPFilteredVoltage

% Derive HP filter (DF) coefficients
% HLP(z) = [1-z^-32] / [1-z^-1]
% H(z) = z^-16 - 1/32*HLP(z)
% H(z) = z^-16 - 1/32 [1-z^-32] / [1-z^-1]
% Y(z)[1-z^-1] = X(z)[1-z^-1]z^-16 - (1/32)X(z)[1-z^-32]
% Y(z)[1-z^-1] = X(z)[z^-16 - z^-17] - X(z)[(1/32) - (1/32)z^-32]
% y(n) - y(n-1) = x(n-16) - x(n-17) - (1/32)x(n) + (1/32)x(n-32)
% Apply to form:
% a(1)y(n)+a(2)y(n-1)=b(1)x(n)+b(2)x(n-1)+b(17)x(n-16)+b(18)x(n-17)+b(33)x(n-32)

HPa = [1, -1];
HPb = zeros(1,33); 
HPb(1) = -(1/32); HPb(17) = 1; HPb(18) = -1; HPb(33) = (1/32);

% Filter the ECG signal
ECGHPFilteredVoltage = filter(HPb, HPa, ECGLPFilteredVoltage);

% ----------------------------------------------------------------- Stage 3
% This stage impliments a Derivative filter
% ECGHPFilteredVoltage -> Derivative filter -> ECGDFFilteredVoltage

% Derive DF coefficients
% H(z) = (1/8) [2 + z^-1 - z^-3 -2z^-4]
% Y(z) = X(z) [ (2/8) + (1/8)z^-1 + (-1/8)z^-3 + (-2/8)z^-4 ]
% y(n) = (2/8)x(n) + (1/8)x(n-1) + (-1/8)x(n-3) + (-2/8)x(n-4)
% Apply to form:
% a(1)y(n)=b(1)x(n)+b(2)x(n-1)+b(3)x(n-2) ect ...

DFa = 1;
DFb = zeros(1,33); 
DFb(1) = (2/8); DFb(2) = (1/8); DFb(4) = (-1/8); DFb(5) = (-2/8);

% Set the first sample as the baseline of the signal
ECGHPFilteredVoltage = ECGHPFilteredVoltage - ECGHPFilteredVoltage(1);

% Filter the ECG signal
ECGDFFilteredVoltage = filter(HPb, HPa, ECGHPFilteredVoltage);

% ----------------------------------------------------------------- Stage 4
% This stage impliments a Squaring Operation
% Input: ECGDFFilteredVoltage
% Output: ECGSquaredVoltage

ECGSquaredVoltage = ECGDFFilteredVoltage.^2;

% ----------------------------------------------------------------- Stage 5
% This stage impliments a 30 Point Moving Average Filter
% Input: ECGSquaredVoltage
% Output: ECGMAVoltage

% Derive 30 pt. MA Filter coefficients
% H(z) = (1/3) [1 + z^-1 + z^-2 + ... + z^-29]
% Y(z) = X(z) (1/3)[1 + z^-1 + z^-2 + ... + z^-29]
% y(n) = (1/3)[x(n) + x(n-1) + x(n-2) + ... + x(n-29)]
% Apply to form:
% a(1)y(n)=b(1)x(n)+b(2)x(n-1)+b(3)x(n-2) ... b(30)x(n-29)
MAa = 1;
MAb = ones(1,30); MAa = MAa / 30;

% Filter the ECG signal
ECGMAVoltage = filter(MAb, MAa, ECGSquaredVoltage);

% ----------------------------------------------------------------- Stage 6
% Blanking
% The code searches through Signal for samples that are larger then their
% local (+35 / -35) neighboring samples.  It stores the location of these
% peaks in PeakIndex.
% Input: ECGMAVoltage = Signal 
% Output: PeakIndex

Signal = ECGMAVoltage; % Call the input signal Signal

Peaks = []; PeakIndexCounter = 0; % Peak is a sample index marking samples that are larger then their local neighbors.

% Loop though the samples in the signal
for GlobalInspectionIndex = 36 : (length(Signal) - 35);
    CurrentSample = Signal(GlobalInspectionIndex);
        
    LocalInspectionIndex = (GlobalInspectionIndex - 35) : (GlobalInspectionIndex + 35);
    MaxLocalSample = max(Signal(LocalInspectionIndex)); % Search through the local samples for a max value
        
    if CurrentSample >= MaxLocalSample % If the current sample is greater then the max local sample
        PeakIndexCounter = PeakIndexCounter + 1;
        PeakIndex(PeakIndexCounter) = GlobalInspectionIndex; % Store the location of the global sample into peaks[];
    end       
    
end

%Plot the the ECG with all peaks marked by red *
PlotPeaks = figure('Name','Peaks Detected - Filtered ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, Signal); 
hold on; plot(ECGTimeAxis(PeakIndex), Signal(PeakIndex), 'r*'); hold off;
title('Peaks Detected - Filtered ECG signal (w/ peaks indicated by red "*") ');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis auto;
subplot(2,1,2); plot(ECGTimeAxis, ECGVoltage); 
hold on; plot(ECGTimeAxis(PeakIndex), ECGVoltage(PeakIndex), 'r*'); hold off;
title('Peaks Detected - Unfiltered ECG signal (w/ peaks indicated by red "*" ');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis auto;

% ------------------------------------------------------------------ Step 7
% Thresholding - This code looks at the peaks of the signal and decides if 
% they are true QRS peaks by checking that they exceed a threshold.  If they 
% are it stores their location in QRSIndex.
% INPUT: Signal & PeakIndex 
% OUTPUT: QRSIndex

% Setup a initial threshold for what we will consider QRS peaks
SignalPower = 280; NoisePower = 57;
Threshold = NoisePower + 0.25 * (SignalPower - NoisePower);

QRSIndex = []; QRSIndexCounter = 0; % An index for QRS Peaks

for PeakIndexCounter = 1: length(PeakIndex) % Loop through all the values of PeakIndex    
    if Signal(PeakIndex(PeakIndexCounter)) > Threshold % If the current sample exceeds the threshold value
                
        SignalPower = (0.125 * Signal(PeakIndex(PeakIndexCounter))) + (0.875 * SignalPower); % Re-weight Signal Power with the signal value > Threshold
         
        QRSIndexCounter = QRSIndexCounter + 1; % Move to the next location of the QRS Index        
        QRSIndex(QRSIndexCounter) = PeakIndex(PeakIndexCounter); % Save the location of this QRS Peak
    else % The current sample did not exceed the threshold value
        NoisePower = (0.125 * Signal(PeakIndex(PeakIndexCounter))) + (0.875 * NoisePower);  % Re-weight the Noise Power with the value < Threshold
    end
    
    % Update the threshold WRT the updated Signal/Noise Power
    Threshold = NoisePower + 0.25 * (SignalPower - NoisePower);   
end

% ----------------------------------------------------------------- Stage 8
% Correct for Non Stationary Signal - Delays are introduced by each filter,
% which in stationary signal we could just subtract for the index and line
% up with the QRS peaks.  However, the ECG is not a stationary signal so so
% this stage searches before and after the MA peak location on the orginal
% ECG signal for a maximum which is the true peak.
% INPUT: QRSIndex
% OUTPUT: QRSIndex

%Plot the ECG with uncorrected QRS peak locations
QRSNoDelayPlot = figure('Name','Corrected and uncorrected QRS Peaks - Unfiltered ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage); 
hold on; plot(ECGTimeAxis(QRSIndex), ECGVoltage(QRSIndex), 'r*'); hold off;
title('The uncorrected QRS peak locations - Unfiltered ECG signal (w/ peaks indicated by red "*" ');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis auto;

% Remove the filter delay by removing 38 samples
QRSIndex = QRSIndex - 38;

% if the first peak starts at less then 30 samples from the start search function will
% try to search a negative array index when it goes back 10, so we will
% throw this value away.
if QRSIndex(1) <= 30
    QRSIndex = QRSIndex(2:end);
end

% if the last peak starts at less then 30 samples from the end the search function will
% try to search a negative array index when it goes back 10, so we will
% throw this value away.
if QRSIndex(length(QRSIndex)) >= (length(ECGVoltage) - 30)
    QRSIndex = QRSIndex(1:end - 1);
end

% Correct for the non-statinarity of the signal by searching neighboring
% samples for a new maximum value
for QRSIndexCounter = 1: length(QRSIndex) % Go though all the values in

    QRSMaxLocation = QRSIndex(QRSIndexCounter); % the location of the uncorrected QRS peak
    CurrentQRSMaxValue = ECGVoltage(QRSIndex(QRSIndexCounter)); % the value at the uncorrected QRS peak
    
    % Search 10 samples on both sides of the uncorrected peak location for
    % the largest value; this is the true peak
    for ECGVoltageLocalCounter = QRSIndex(QRSIndexCounter) - 30 : QRSIndex(QRSIndexCounter) + 30
                   
        if CurrentQRSMaxValue < ECGVoltage(ECGVoltageLocalCounter); % If we find a new maximum value            
            CurrentQRSMaxValue = ECGVoltage(ECGVoltageLocalCounter); % Save the new max value
            QRSMaxLocation = ECGVoltageLocalCounter; % Save the new max location
        end    
    
    end
    
    QRSIndex(QRSIndexCounter) = QRSMaxLocation; % Update the QRS Peak location index with the new max location    
 end

% Plot the corrected QRS peak locations
subplot(2,1,2); plot(ECGTimeAxis, ECGVoltage);
hold on; plot(ECGTimeAxis(QRSIndex), ECGVoltage(QRSIndex), 'r*'); hold off;
title('The corrected QRS peak locations - Unfiltered ECG signal (w/ peaks indicated by red "*" ');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis auto;