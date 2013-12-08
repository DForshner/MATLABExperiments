% Filtering of the ECG for the Removal of Noise
%
% The algorithm developed by Pan and Tompkins identiﬁes QRS complexes based 
% on analysis of the slope, amplitude, and width of the QRS. 

close all; clear all; clc; % Clear everything
pause on; % enable pausing

% --------------------------------------------------------------- Load File
load ecgpvcs; % ecgpvcs has two signals: ecg1 and ecg2
T = 10 / 2000 % 2,000 samples = 10 seconds 
fs = 1 / T

% ------------------------------------------------------------------------
% ----------------------------------------------- Part 1: Training Section
% ------------------------------------------------------------------------

% ------------------------------------------------------------------ Step 1
%1. Use the ﬁrst 40% of the given ECG signal for the training step.
TrainingECG = ecg1(1:round( 0.40 * length(ecg1)));
TrainingECGTimeAxis = ( 1:length(TrainingECG) ) * T;

% ------------------------------------------------------------------ Step 2
%2. Apply the Pan–Tompkins method and detect the beats in the signal.
QRSPeakIndex = FindQRS(TrainingECG, TrainingECGTimeAxis, fs);

% ------------------------------------------------------------------ Step 3
%3. Segment each cardiac cycle (QRS-T complex) by taking a few samples before and
%   a few samples after the corresponding marker point detected in the output of the
%   Pan–Tompkins method. Use a total duration of about 300 ms for each beat. The
%   P wave need not be considered in the present exercise. Prepare a plot of the entire
%   ECG signal with the segmentation points marked for each beat in the signal.
SegLength = round(0.300 / T)

for QRSPeakIndexCounter = 1:length(QRSPeakIndex)        
    
    StartSegmentIndex(QRSPeakIndexCounter) = QRSPeakIndex(QRSPeakIndexCounter) - 10;    
    EndSegmentIndex(QRSPeakIndexCounter) = StartSegmentIndex(QRSPeakIndexCounter) + SegLength;
    
 end

TrainECGSegmented = figure('Name','Training: ECG with segments marked - Unfiltered ECG signal'); % Create a new figure
plot(TrainingECGTimeAxis, TrainingECG); 
hold on; plot(TrainingECGTimeAxis(StartSegmentIndex), TrainingECG(StartSegmentIndex), 'r*');
plot(TrainingECGTimeAxis(EndSegmentIndex), TrainingECG(EndSegmentIndex), 'r+'); hold off;
title('Training: ECG with segments marked - Unfiltered ECG signal (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis auto;

% ------------------------------------------------------------------ Step 4
%4. Prepare an array with one line per beat in the training signal, giving the ECG beat
%   number; the type of the beat as normal (0) or PVC (1); QRST A; and CC . Save this
%   array for use in the testing step.
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    TrainingArray(QRSPeakIndexCounter, 1) = QRSPeakIndexCounter; % ECG Beat number
    
    TrainingArray(QRSPeakIndexCounter, 2) = FindRRInterval(QRSPeakIndex, QRSPeakIndexCounter, fs); % RR Interval
        
    TrainingArray(QRSPeakIndexCounter, 3) = FindQRSTArea(TrainingECG, TrainingECGTimeAxis, fs, StartSegmentIndex(QRSPeakIndexCounter), EndSegmentIndex(QRSPeakIndexCounter)); 
                
    Template = TrainingECG(StartSegmentIndex(1) : EndSegmentIndex(1) );
    Signal = TrainingECG(StartSegmentIndex(QRSPeakIndexCounter) : EndSegmentIndex(QRSPeakIndexCounter));
    TrainingArray(QRSPeakIndexCounter, 4) = FindCC(Template, Signal); % Cross Corralation               
    
    % Decision line points (0, 25) to (1, 100) so y = x(75) + 25    
    if TrainingArray(QRSPeakIndexCounter, 3) >= (TrainingArray(QRSPeakIndexCounter, 2) * 75 + 25)
        TrainingArray(QRSPeakIndexCounter, 5) = 1; % PVC
    else
        TrainingArray(QRSPeakIndexCounter, 5) = 0; % Normal
    end
       
end

% The ranges of the values of CC and QRST A could be substantially diﬀerent,
% which creates diﬃculties in k −NN analysis. To overcome this, normalize the QRST A
% values in the training set and also in the testing set by dividing by the maximum value of
% the feature in the training set.
TrainingArray(:, 3) = TrainingArray(:, 3) / max(TrainingArray(:, 3));

% ------------------------------------------------------------------ Step 5
%5. Prepare a plot of the training portion of the given signal, with each beat labeled as
%   ‘o’ or ‘x’, representing a normal beat or PVC, respectively.

TrainECGTyped = figure('Name','Training: ECG with beats types indicated - Unfiltered ECG signal'); % Create a new figure
subplot(4,1,1); plot(TrainingECGTimeAxis, TrainingECG); hold on;
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)   
    if TrainingArray(QRSPeakIndexCounter, 5) == 0
        plot(TrainingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TrainingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rO');
    else
        plot(TrainingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TrainingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rX');
    end    
end
hold off;
title('Training: ECG with beats types indicated - Unfiltered ECG signal (w/ Normal "o" and PVC "x")');
ylabel('ECG (mV)'); axis ([0 25 1400 3500]);

subplot(4,1,2); plot(TrainingECGTimeAxis, TrainingECG); hold on;
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)   
    if TrainingArray(QRSPeakIndexCounter, 5) == 0
        plot(TrainingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TrainingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rO');
    else
        plot(TrainingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TrainingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rX');
    end    
end
hold off;
ylabel('ECG (mV)'); axis ([25 50 1400 3500]);

subplot(4,1,3); plot(TrainingECGTimeAxis, TrainingECG); hold on;
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)   
    if TrainingArray(QRSPeakIndexCounter, 5) == 0
        plot(TrainingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TrainingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rO');
    else
        plot(TrainingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TrainingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rX');
    end    
end
hold off;
ylabel('ECG (mV)'); axis ([50 75 1400 3500]);

subplot(4,1,4); plot(TrainingECGTimeAxis, TrainingECG); hold on;
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)   
    if TrainingArray(QRSPeakIndexCounter, 5) == 0
        plot(TrainingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TrainingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rO');
    else
        plot(TrainingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TrainingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rX');
    end    
end
hold off;
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis ([75 100 1400 3500]);

% ------------------------------------------------------------------ Step 6
%6. Prepare a scatter plot of CC versus QRST A for all of the beats detected. Mark the
%   plot with ‘o’ or ‘x’ for each normal beat or PVC, respectively.
TrainScatterPlot = figure('Name','Scatter plot of CC vs. QRST A for all of the beats detected'); % Create a new figure
scatter(0, 0, 'wo'); hold on; % mark origin with white dot

for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    if TrainingArray(QRSPeakIndexCounter, 5) == 0
        scatter(TrainingArray(QRSPeakIndexCounter,2), TrainingArray(QRSPeakIndexCounter,3), 'bo');
    else
        scatter(TrainingArray(QRSPeakIndexCounter,2), TrainingArray(QRSPeakIndexCounter,3), 'rx');
    end
    
end
hold off;
title('Training: CC vs. QRST A for all of the beats detected (w/ Normal "o" and PVC "x")');
xlabel('CC'); ylabel('QRST Area'); axis auto;


% ------------------------------------------------------------------------
% ------------------------------------------------ Part 2: Testing Section
% ------------------------------------------------------------------------

% ------------------------------------------------------------------ Step 1
%1. Apply the QRS detection and parameterization procedures as above to the remaining
%   part of the ECG signal.

TestingECG = ecg1;
TestingECGTimeAxis = (1:length(ecg1)) * T;

N = round(0.4 * length(TestingECG));
TestingECG = TestingECG( N : end);
TestingECGTimeAxis = TestingECGTimeAxis( N : end);

QRSPeakIndex = FindQRS(TestingECG, TestingECGTimeAxis, fs);

SegLength = round(0.300 / T)

for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    StartSegmentIndex(QRSPeakIndexCounter) = QRSPeakIndex(QRSPeakIndexCounter) - 10;    
    EndSegmentIndex(QRSPeakIndexCounter) = StartSegmentIndex(QRSPeakIndexCounter) + SegLength;
    
 end

TestECGSegmented = figure('Name','Testing: ECG with segments marked - Unfiltered ECG signal'); % Create a new figure
plot(TestingECGTimeAxis, TestingECG); 
hold on; plot(TestingECGTimeAxis(StartSegmentIndex), TestingECG(StartSegmentIndex), 'r*');
plot(TestingECGTimeAxis(EndSegmentIndex), TestingECG(EndSegmentIndex), 'r+'); hold off;
title('Testing: ECG with segments marked - Unfiltered ECG signal (w/ segment starts red "*" ends red "+"');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis auto;

% ------------------------------------------------------------------ Step 2
%2. For each beat detected, form a feature vector as [CC, QRST A].
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    TestingArray(QRSPeakIndexCounter, 1) = QRSPeakIndexCounter; % ECG Beat number
    
    TestingArray(QRSPeakIndexCounter, 2) = FindRRInterval(QRSPeakIndex, QRSPeakIndexCounter, fs); % RR Interval
        
    TestingArray(QRSPeakIndexCounter, 3) = FindQRSTArea(TestingECG, TestingECGTimeAxis, fs, StartSegmentIndex(QRSPeakIndexCounter), EndSegmentIndex(QRSPeakIndexCounter)); 
                
    Template = TestingECG(StartSegmentIndex(1) : EndSegmentIndex(1) );
    Signal = TestingECG(StartSegmentIndex(QRSPeakIndexCounter) : EndSegmentIndex(QRSPeakIndexCounter));
    TestingArray(QRSPeakIndexCounter, 4) = FindCC(Template, Signal); % Cross Corralation               
    
    TestingArray(QRSPeakIndexCounter, 5) = 3; % 3 is error will be defined later by k nearest neighbor    
       
end

% ------------------------------------------------------------------ Step 3
%3. For each beat, ﬁnd the k nearest neighbors in the training set using the feature vector
%   and the Euclidean distance measure, and classify the beat as normal or PVC using
%   the k −NN rule with k = 1. (See Section 9.4.3 of the textbook.)
%   Repeat the procedure with k = 3.

% The ranges of the values of CC and QRST A could be substantially diﬀerent,
% which creates diﬃculties in k −NN analysis. To overcome this, normalize the QRST A
% values in the training set and also in the testing set by dividing by the maximum value of
% the feature in the training set.
TestingArray(:, 3) = TestingArray(:, 3) / max(TestingArray(:, 3));

% ----------------------------------------------------------- K = 1
% Then for each beat in the testing set find the nearest neighbor (for K=1) 
% and 3 nearest neighbors (k=3) (that is, calculate
% ED for each beat in the testing set to all beats in the training set, sort
% them, and find the nearest neighbor type). If the nearest neighbor is
% abnormal, you classify that beat as abnormal. For K=3, you check if 2
% nearest neighbors are abnormal.  Finally, you manually mark each beat in the
% test set (Figure 9.7 in book), and report the TP, FP, accuracy etc.

EDistance = [];

for TestingArrayCounter = 1:length(TestingArray)
 
    % We are going to examine every point in the training array to find ED
    for TrainingArrayCounter = 1:length(TrainingArray)                
        Term1 = ( TestingArray(TestingArrayCounter, 2) - TrainingArray(TrainingArrayCounter, 2) ).^2;
        Term2 = ( TestingArray(TestingArrayCounter, 3) - TrainingArray(TrainingArrayCounter, 3) ).^2;
        EDistance(TrainingArrayCounter) = sqrt(Term1 + Term2);
    end
    
    % Now locate the closest point and see what type it is    
    CurrentMinED = EDistance(1);
    CurrentMinIndex = 1;
    
    for TrainingArrayCounter = 1:length(TrainingArray)                
        
        if EDistance(TrainingArrayCounter) <= CurrentMinED            
            CurrentMinED = EDistance(TrainingArrayCounter);
            CurrentMinIndex = TrainingArrayCounter;
        end    
                
    end
    
    % The test sample is the same as it's closest neighbor
    TestingArray(TestingArrayCounter, 5) = TrainingArray(CurrentMinIndex, 5);    
    
end

TestingScatterPlotK1 = figure('Name','Testing (K=1) : Scatter plot of CC vs. QRST A for all of the beats detected'); % Create a new figure
scatter(0, 0, 'wo'); hold on; % mark origin with white dot

for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    if TestingArray(QRSPeakIndexCounter, 5) == 0
        scatter(TestingArray(QRSPeakIndexCounter,2), TestingArray(QRSPeakIndexCounter,3), 'bo');
    else
        scatter(TestingArray(QRSPeakIndexCounter,2), TestingArray(QRSPeakIndexCounter,3), 'rx');
    end
    
end
hold off;
title('Testing: CC vs. QRST A for all of the beats detected (w/ Normal "o" and PVC "x")');
xlabel('CC'); ylabel('QRST Area'); axis auto;

% ----------------------------------------------------------- K = 3
EDistance = [];

for TestingArrayCounter = 1:length(TestingArray)
 
    % We are going to examine every point in the training array to find ED
    for TrainingArrayCounter = 1:length(TrainingArray)                
        Term1 = ( TestingArray(TestingArrayCounter, 2) - TrainingArray(TrainingArrayCounter, 2) ).^2;
        Term2 = ( TestingArray(TestingArrayCounter, 3) - TrainingArray(TrainingArrayCounter, 3) ).^2;
        EDistance(TrainingArrayCounter) = sqrt(Term1 + Term2);
    end
    
    % Now locate the 3 closest points and see what types they are
       
    [SortED,SortEDIndex] = sort(EDistance);
       
    NearestThreeSum = sum(TrainingArray(SortEDIndex(1), 5) + TrainingArray(SortEDIndex(2), 5) + TrainingArray(SortEDIndex(3), 5))    
       
    if NearestThreeSum == 3
       TestingArray(TestingArrayCounter, 5) = 1; % Is PVC    
    elseif NearestThreeSum == 2
       TestingArray(TestingArrayCounter, 5) = 1; % Is PVC
    else
       TestingArray(TestingArrayCounter, 5) = 0; % Is Normal
    end    
           
end

TestingScatterPlotK3 = figure('Name','Testing (K=3) : Scatter plot of CC vs. QRST A for all of the beats detected'); % Create a new figure
scatter(0, 0, 'wo'); hold on; % mark origin with white dot

for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    if TestingArray(QRSPeakIndexCounter, 5) == 0
        scatter(TestingArray(QRSPeakIndexCounter,2), TestingArray(QRSPeakIndexCounter,3), 'bo');
    else
        scatter(TestingArray(QRSPeakIndexCounter,2), TestingArray(QRSPeakIndexCounter,3), 'rx');
    end
    
end
hold off;
title('Testing: CC vs. QRST A for all of the beats detected (w/ Normal "o" and PVC "x")');
xlabel('CC'); ylabel('QRST Area'); axis auto;

% ------------------------------------------------------------------ Step 4
%4. Prepare a plot of the testing portion of the given signal, with each beat labeled as ‘o’
%   or ‘x’, representing a normal beat or PVC, respectively, as determined by the k −NN
%   rule in your program.

TestECGTyped = figure('Name','Testing: ECG with beats types indicated - Unfiltered ECG signal'); % Create a new figure
subplot(4,1,1); plot(TestingECGTimeAxis, TestingECG); hold on;
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    if TestingArray(QRSPeakIndexCounter, 5) == 0
        plot(TestingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TestingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rO');
    else
        plot(TestingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TestingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rX');
    end
    
end
hold off; title('Testing: ECG with beats types indicated - Unfiltered ECG signal (w/ Normal "o" and PVC "x")');
ylabel('ECG (mV)'); axis ([100 138 1400 3700]);

subplot(4,1,2); plot(TestingECGTimeAxis, TestingECG); hold on;
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    if TestingArray(QRSPeakIndexCounter, 5) == 0
        plot(TestingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TestingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rO');
    else
        plot(TestingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TestingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rX');
    end
    
end
hold off; title('Testing: ECG with beats types indicated - Unfiltered ECG signal (w/ Normal "o" and PVC "x")');
ylabel('ECG (mV)'); axis ([138 176 1400 3700]);

subplot(4,1,3); plot(TestingECGTimeAxis, TestingECG); hold on;
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    if TestingArray(QRSPeakIndexCounter, 5) == 0
        plot(TestingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TestingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rO');
    else
        plot(TestingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TestingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rX');
    end
    
end
hold off; title('Testing: ECG with beats types indicated - Unfiltered ECG signal (w/ Normal "o" and PVC "x")');
ylabel('ECG (mV)'); axis ([176 214 1400 3700]);

subplot(4,1,4); plot(TestingECGTimeAxis, TestingECG); hold on;
for QRSPeakIndexCounter = 1:length(QRSPeakIndex)
    
    if TestingArray(QRSPeakIndexCounter, 5) == 0
        plot(TestingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TestingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rO');
    else
        plot(TestingECGTimeAxis(QRSPeakIndex(QRSPeakIndexCounter)), TestingECG(QRSPeakIndex(QRSPeakIndexCounter)), 'rX');
    end
    
end
hold off; title('Testing: ECG with beats types indicated - Unfiltered ECG signal (w/ Normal "o" and PVC "x")');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis ([214 250 1400 3700]);

% ------------------------------------------------------------------ Step 5
%5. Check the results of your classiﬁcation procedure, and compute the accuracy of clas-
%   siﬁcation, including measures of true-positive, false-positive, true-negative, and false-
%   negative fractions. Note also the number of beats not detected (if any) by your
%   program. Prepare a summary table of your results.



% ------------------------------------------------------------------
% Calcaulate the Mean and StdDev of CC & QRSTArea

NormalBeatCounter = 1;
PVCBeatCounter = 1;
for TrainingArrayIndexCounter = 1:length(TrainingArray)
    
    if TrainingArray(TrainingArrayIndexCounter, 5) == 0        
        TrainNormalBeatCC(NormalBeatCounter) = TrainingArray(TrainingArrayIndexCounter, 4);
        TrainNormalBeatQRSTA(NormalBeatCounter) = TrainingArray(TrainingArrayIndexCounter, 3);
        NormalBeatCounter = NormalBeatCounter + 1;
    else
        TrainPVCBeatCC(PVCBeatCounter) = TrainingArray(TrainingArrayIndexCounter, 4);
        TrainPVCBeatQRSTA(PVCBeatCounter) = TrainingArray(TrainingArrayIndexCounter, 3);
        PVCBeatCounter = PVCBeatCounter + 1;
    end
    
end
TrainNormalBeatCCmean = mean(TrainNormalBeatCC)
TrainNormalBeatCCstd = std(TrainNormalBeatCC)
TrainNormalBeatQRSTAmean = mean(TrainNormalBeatQRSTA)
TrainNormalBeatQRSTAstd = std(TrainNormalBeatQRSTA)

mean(TrainPVCBeatCC)
std(TrainPVCBeatCC)
mean(TrainPVCBeatQRSTA)
std(TrainPVCBeatQRSTA)

NormalBeatCounter = 1;
PVCBeatCounter = 1;
for TestingArrayIndexCounter = 1:length(TestingArray)
    
    if TestingArray(TestingArrayIndexCounter, 5) == 0        
        NormalBeatCC(NormalBeatCounter) = TestingArray(TestingArrayIndexCounter, 4);
        NormalBeatQRSTA(NormalBeatCounter) = TestingArray(TestingArrayIndexCounter, 3);
        NormalBeatCounter = NormalBeatCounter + 1;
    else
        PVCBeatCC(PVCBeatCounter) = TestingArray(TestingArrayIndexCounter, 4);
        PVCBeatQRSTA(PVCBeatCounter) = TestingArray(TestingArrayIndexCounter, 3);
        PVCBeatCounter = PVCBeatCounter + 1;
    end
    
end
mean(NormalBeatCC)
std(NormalBeatCC)
mean(NormalBeatQRSTA)
std(NormalBeatQRSTA)

mean(PVCBeatCC)
std(PVCBeatCC)
mean(PVCBeatQRSTA)
std(PVCBeatQRSTA)
