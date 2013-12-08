% Analysis of the Relationship Between Parameters of the EMG Signal and Muscular Force.
% Using the EMG data we recorded in lab# 1, this program will calculate a
% series of parameters for the EMG signal and muscular force.

close all; % Clear all windows
clear all; % Clear all variables
clc; % Clear console

load EMGAndForceVsTime.txt;  % Load file into memory
EMG_Time = EMGAndForceVsTime(:,1);   % copy column 1 of the matrix to Time
EMG_Force = EMGAndForceVsTime(:,2);   % copy column 2 of the matrix to Voltage
EMG_Voltage = EMGAndForceVsTime(:,3);   % copy column 3 of the matrix to Force

% ---------------------------------------------------------- Part 1
% Normalize the Force vector
EMG_ForceShifted = EMG_Force - min(EMG_Force);
EMG_ForceNorm = EMG_ForceShifted / max(EMG_ForceShifted) * 100; % Note that min(EMG_ForceShifted is zero)

%Subtract the mean from the EMG signal to remove any bias
EMG_Voltage = EMG_Voltage - mean(EMG_Voltage);

% ---------------------------------------------------------- Part 2
% Manually identify portions (segments) corresponding to each level of
% contraction within which the force remains almost constant
IdentSegTime = [4.0, 6.2, 9.8, 12.6, 17.1, 20.1, 24.5, 28.7, 33.9, 37.8, 42.5, 47, 52.6 ,57.1];

% Plot of the EMG (mV) and normalized force (in percent % MVC) vs  Time
% (Sec) with segments marked
EMG_VoltageForceVsTimeSeg=figure('Name','EMG signal and Normalized force vs. Time'); % Create a new figure
subplot(2,1,1); plot(EMG_Time, EMG_ForceNorm)
hold on;
plot(IdentSegTime,0,'r*');  % Highlight 
hold off;
ylabel('Normalized Force (%MVC)');
xlabel('Time (Sec)');
subplot(2,1,2); plot(EMG_Time, EMG_Voltage)
hold on;
plot(IdentSegTime,0,'r*');  % Highlight 
hold off;
ylabel('EMG Signal (mV)');
xlabel('Time (Sec)');
title('EMG signal and Normalized force vs. Time');
axis auto; 

% ---------------------------------------------------------- Part 3
% For each segment of the EMG identified compute the DR, MS , RMS and ZCR
% parameters.

SampleInterval = 0.0005; % The EMG signals are sampled with a period of 0.0005s

%ForceSegment = figure; % uncomment to see force intervals

for Counter=1:1:7
    Segment = Counter * 2 - 1;
    disp(['Segment = ', num2str(Segment)]); % Display in console
    disp(['Start = ', num2str(IdentSegTime(Segment))]); % Display in console
    disp(['Stop = ', num2str(IdentSegTime(Segment + 1))]); % Display in console
    
    % Identify the segments of interest from the EMG signals
    s1 = IdentSegTime(Segment) / SampleInterval; % First point index
    s2 = IdentSegTime(Segment + 1) / SampleInterval; % Second point index
    EMG_ForceSeg = EMG_ForceNorm(s1:s2); % Segment of interest of force signal
    EMG_VoltageSeg = EMG_Voltage(s1:s2); % Segment of interest of voltage signal  
        
    % Calculate Dynamic Range
    DynamicRange(Counter) = max(EMG_VoltageSeg) - min(EMG_VoltageSeg);
       
    % Calculate Mean Squared Value
    NumSamplesInSeg = s2 - s1 + 1; % The number of samples in the segment
    MeanSquared(Counter) = sum(EMG_VoltageSeg.^2) / NumSamplesInSeg;
        
    % Calculate RMS
    RMS(Counter) = sqrt(MeanSquared(Counter));
    
    % Calculate ZCR
    NumZC = 0; % Number of zero crossings    
    for ZCRcounter = 1:(s2 - s1)       
       if (sign(EMG_VoltageSeg(ZCRcounter)) + sign(EMG_VoltageSeg(ZCRcounter + 1))) == 0
            NumZC = NumZC + 1;
       end
    end
    ZCR(Counter) = NumZC / ((s2 - s1 + 1) * SampleInterval);
    
    % Calculate mean voltage value for each segment
    ForceMean(Counter) = mean(EMG_ForceSeg)
    
    %subplot(14,1,Segment); plot(EMG_ForceSeg); %uncomment plots force segments        
end

% ---------------------------------------------------------- Part 4
% Plot the DR, MS, RMS, and ZCR values versus force in %MVC.  Lable the
% axis with appropriate units.

DRMSRMSZCRvsMVC=figure('Name','EMG parameters vs. force'); % Create a new figure
subplot(4,1,1); plot(ForceMean, DynamicRange)
title('Dynamic Range vs. Normalized Force');
ylabel('DR (mV)');
xlabel('Normalized Force (%MVC)');
subplot(4,1,2); plot(ForceMean, MeanSquared)
title('Mean Squared (Average Power) vs. Normalized Force');
ylabel('MS (mV^2)');
xlabel('Normalized Force (%MVC)');
subplot(4,1,3); plot(ForceMean, RMS)
title('Root Mean Squared (average magnitude) vs. Normalized Force');
ylabel('RMS (mV)');
xlabel('Normalized Force (%MVC)');
subplot(4,1,4); plot(ForceMean, ZCR)
title('Zero Crossing Rate vs. Normalized Force');
ylabel('ZCR (# / sec)');
xlabel('Normalized Force (%MVC)');
axis auto; 

% ---------------------------------------------------------- Part 5
% Using the polyfit function obtain a straight line fit to represent teh
% variation of each EMG parametere vs. Force.  Use polyval to evaluate the
% values of the dependant variable giving by the models and plot them vs.
% the orignal values

Counter = (1:1:7); % segment index variable for linear parameters
% Remember that ForceMean is a size 7 vector of mean force values per seg

% Calculate a linear model for Dynamic Range and add it to the plot
PDynamicRange = polyfit(ForceMean, DynamicRange,1)
subplot(4,1,1); hold on; plot(ForceMean, POLYVAL(PDynamicRange,ForceMean),'r--'); hold off;
title('Dynamic Range (blue) & Linear Model (red) vs. Normalized Force');

% Calculate a linear model for mean squared and add it to the plot
PMeanSquared = polyfit(ForceMean, MeanSquared, 1)
subplot(4,1,2); hold on; plot(ForceMean, POLYVAL(PMeanSquared,ForceMean),'r--'); hold off;
title('Mean Squared (Average Power) (blue) & Linear Model (red) vs. Normalized Force');

% Calculate a linear model for Root mean squared and add it to the plot
PRMS = polyfit(ForceMean, RMS, 1)
subplot(4,1,3); hold on; plot(ForceMean, POLYVAL(PRMS,ForceMean),'r--'); hold off;
title('Root Mean Squared (average magnitude) (blue) & Linear Model (red) vs. Normalized Force');

% Calculate a linear model for Zero Crossing rate and add it to the plot
PZCR = polyfit(ForceMean, ZCR, 1)
subplot(4,1,4); hold on; plot(ForceMean, POLYVAL(PZCR,ForceMean),'r--'); hold off;
title('Zero Crossing Rate (blue) & Linear Model (red) vs. Normalized Force');

% ---------------------------------------------------------- Part 6
% Computer the correlation coefficient r^2 and analyze the goodness of fit
% for each parameter

% Calculate the correlation coefficient for DR
x = ForceMean; % calculated parameter
%y = POLYVAL(PDynamicRange, ForceMean) % Linear Model
y = DynamicRange;
N = 7; % Index size

%Verify formula to calculate the top works
%Top = 0;
%for i = 1:N 
%    Top = Top + ( x(i) * y(i) )
%end
%Top = (Top - N * mean(x) * mean(y) ).^2
Top = ( sum(x .* y) - N * mean(x) * mean(y) ).^2

%Verify formula to calculate the bottomL works
%BottomL = 0;
%for i = 1:N 
%    BottomL = BottomL + ( x(i).^2 )    
%end
%BottomL = BottomL - N * mean(x).^2

BottomL = sum(x.^2) - N * ((mean(x)).^2)
BottomR = sum(y.^2) - N * ((mean(y)).^2)
RRDynamicRange = Top / (BottomL * BottomR)

% Calculate the correlation coefficient for MS
x = ForceMean; % calculated parameter
%y = POLYVAL(PMeanSquared, ForceMean); % Linear Model
y = MeanSquared;
N = 7; % Index size
Top = ( sum(x .* y) - N * mean(x) * mean(y) ).^2;
BottomL = sum(x.^2) - N * (mean(x).^2);
BottomR = sum(y.^2) - N * (mean(y).^2);
RRMeanSquared = Top / (BottomL * BottomR)

% Calculate the correlation coefficient for RMS
x = ForceMean; % calculated parameter
%y = POLYVAL(PRMS, ForceMean) % Linear Model
y = RMS;
N = 7; % Index size
Top = ( sum(x .* y) - N * mean(x) * mean(y) ).^2;
BottomL = sum(x.^2) - N * (mean(x).^2);
BottomR = sum(y.^2) - N * (mean(y).^2);
RRRMS = Top / (BottomL * BottomR)

% Calculate the correlation coefficient for ZCR
x = ForceMean; % calculated parameter
%y = POLYVAL(PZCR, ForceMean) % Linear Model
y = ZCR;
N = 7; % Index size
Top = ( sum(x .* y) - N * mean(x) * mean(y) ).^2;
BottomL = sum(x.^2) - N * (mean(x).^2);
BottomR = sum(y.^2) - N * (mean(y).^2);
RRZRC = Top / (BottomL * BottomR)

% ---------------------------------------------------------- Results
% Tabulate the measured and calculated data

% Tabulate the measured parameters and mean force
ParameterTable = figure('Name','Parameters of the linear model'); % Create a new figure
data = [ForceMean' POLYVAL(PDynamicRange, ForceMean)' POLYVAL(PMeanSquared, ForceMean)' POLYVAL(PRMS, ForceMean)' POLYVAL(PZCR, ForceMean)'];
colnames = {'Mean Force' 'Dynamic Range', 'Mean Squared', 'Root Mean Squared', 'Zero Crossing Rate'};
rownames = {'segment 1', 'segment 2', 'segment 3', 'segment 4', 'segment 5', 'segment 6', 'segment 7'};
t = uitable(ParameterTable, 'Data', data, 'ColumnName', colnames, 'RowName', rownames, 'Units', 'Normalized', 'Position',[0 0 1 1]);

% Tabulate the parameters of the linear model and r^2 for eac variable
ParameterTable = figure('Name','Linear model parameters and correlation coefficient'); % Create a new figure
data = [PDynamicRange' PMeanSquared' PRMS' PZCR' ; [RRDynamicRange RRMeanSquared RRRMS RRZRC] ];
colnames = {'Dynamic Range', 'Mean Squared', 'Root Mean Squared', 'Zero Crossing Rate'};
rownames = {'P1', 'P2', 'r^2'};
t = uitable(ParameterTable, 'Data', data, 'ColumnName', colnames, 'RowName', rownames, 'Units', 'Normalized', 'Position',[0 0 1 1]);
