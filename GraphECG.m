% 2.3 Plot your ECG signals
% Loads the ECG data files in vectors and displays 10 seconds worth of each
% shows a zoomed 2 cycles of the clean ECG.

clear all % Clear all Varables
close all % Clear all windows

ThreeECGs=figure('Name','Three types of ECGs from Lab 1'); % Create a new figuure for the 3 ECGs

load CleanECG.txt;  % Load file into memory
% We have a sample Interval=0.005s so I need 2000 sample points for 10s
CleanECG_Time = CleanECG(1:2000,1);   %  copy column 1 of the matrix to Time
CleanECG_Voltage = CleanECG(1:2000,2);   %  copy column 2 of the matrix to Voltage
subplot(3,1,1), plot(CleanECG_Time,CleanECG_Voltage)
title('ECG #1 - clean signal');
xlabel('Time Frame (sec)');
ylabel('Voltage (mv)');
axis([0 10 -1.5 1.5]);

load MessyECG.txt;  % Load file into memory
% We have a sample Interval=0.005s so I need 2000 sample points for 10s
MessyECG_Time = MessyECG(1:2000,1);   %  copy column 1 of the matrix to Time
MessyECG_Voltage = MessyECG(1:2000,2);   %  copy column 2 of the matrix to Voltage
subplot(3,1,2), plot(MessyECG_Time,MessyECG_Voltage)
title('ECG #2 - noisy signal at a certain heart rate');
xlabel('Time Frame (sec)');
ylabel('Voltage (mv)');
axis([0 10 -1.5 1.5]);

load HigherECG.txt;  % Load file into memory
% We have a sample Interval=0.005s so I need 2000 sample points for 10s
HigherECG_Time = HigherECG(2000:4000,1);   %  copy column 1 of the matrix to Time
HigherECG_Voltage = HigherECG(2000:4000,2);   %  copy column 2 of the matrix to Voltage
subplot(3,1,3), plot(HigherECG_Time,HigherECG_Voltage)
title('ECG #3 - noisy signal at a higher heart rate');
xlabel('Time Frame (sec)');
ylabel('Voltage (mv)');
axis([10 20 -1.5 1.5]);

TwoCyclesECG=figure('Name','Two full cardiac cycles of the clean ECG signal'); % Create a new figuure for the 2 Cycle ECG display
TwoCyclesCleanECG_Time = CleanECG(1:400,1);   %  copy column 1 of the matrix to Time
TwoCyclesCleanECG_Voltage = CleanECG(1:400,2);   %  copy column 2 of the matrix to Voltage
plot(TwoCyclesCleanECG_Time,TwoCyclesCleanECG_Voltage)
title('ECG - 2 Complete Cycles of the clean signal');
xlabel('Time Frame (sec)');
ylabel('Voltage (mv)');
axis([0 2 -1 1]);
