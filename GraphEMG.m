% 2.4 Plot your EMG/force signals
% Plots the EMG/force data we acquired in lab 1.  Creates a second plot
% showing a zoomed in portion of the EMG at the beginning of a contaction.

clear all % Clear all Varables
close all % Clear all windows

%  Create a two-channel plot of the EMG/force data acquired in Lab 1

EMGAndForceVsTime=figure('Name','Two-channel EMG and Force vs. Time'); % Create a new figure

load EMGAndForceVsTime.txt;  % Load file into memory
EMGAndForceVsTime_Time = EMGAndForceVsTime(:,1);   % copy column 1 of the matrix to Time
EMGAndForceVsTime_Voltage = EMGAndForceVsTime(:,2);   % copy column 2 of the matrix to Voltage
EMGAndForceVsTime_Force = EMGAndForceVsTime(:,3);   % copy column 3 of the matrix to Force

subplot(2,1,1); plot(EMGAndForceVsTime_Time,EMGAndForceVsTime_Voltage)
title('Force vs. Time');
xlabel('Time Frame (sec)');
ylabel('Voltage (mv)');
axis([0 60 -12 -8.5]);

subplot(2,1,2); plot(EMGAndForceVsTime_Time,EMGAndForceVsTime_Force)
title('EMG vs. Time');
xlabel('Time Frame (sec)');
ylabel('Voltage (mv)');
axis([0 60 -3 3]);

% Create a two-channel plot expanded showing a segment of about 200 %milliseconds 
% with increasing force at the beginning of one of the contractions and the corresponding EMG.

% Create a new figure for the 2 Cycle ECG display
ExpandedEMGAndForceVsTime=figure('Name','Expanded two-channel EMG and Force vs. Time plot');

subplot(2,1,1); plot(EMGAndForceVsTime_Time,EMGAndForceVsTime_Voltage)
title('Force vs. Time');
xlabel('Time Frame (sec)');
ylabel('Voltage (mv)');
axis([52.2 52.4 -12 -8]); % Use axis to single out section of interest

subplot(2,1,2); plot(EMGAndForceVsTime_Time,EMGAndForceVsTime_Force)
title('EMG vs. Time');
xlabel('Time Frame (sec)');
ylabel('Voltage (mv)');
axis([52.2 52.4 -1.65 1.65]);  % Use axis to single out section of interest