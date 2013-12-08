% Filtering of the ECG for the Removal of Noise
%
% Part 1 - Apply the von Hann lowpass ﬁlter.  Specify the
% ﬁlter in terms of the a and b arrays via the filter command in MATLAB.
% Obtain its freqency responce (magnitude and phase), pole-zero plot, 
% as well as the Fourier spectra of the input and output signals.
%
% Part 2 - Modify the derivative-based ﬁlter given by Equation 3.47 in the textbook for the
% removal of low-frequency artifacts so that the gain at the maximum frequency present
% in the input signal is unity. Use the filter command and apply the ﬁlter to your
% signal.  Obtain its freqency responce (magnitude and phase), pole-zero plot, 
% as well as the Fourier spectra of the input and output signals.
% 
% Part 3 - Apply the notch ﬁlter that you designed in Lab 6 for the rejection of 60 Hz to your
% signal. Ensure that the ﬁlter is normalized to have unit gain at DC. Use the filter
% command.  Obtain its freqency responce (magnitude and phase), pole-zero plot, 
% as well as the Fourier spectra of the input and output signals.
%
% Part 4 - Apply all three ﬁlters to the ECG signal in series, and study the combined ﬁlter and
% the result as speciﬁed above.  Obtain its freqency responce (magnitude and phase), 
% pole-zero plot, as well as the Fourier spectra of the input and output signals.
%
% Part 5 - Apply all three filters by finding the impulse response of the combined ﬁlter 
% as the convolution of the impulse responses of the individual ﬁlters.

close all; clear all; clc; % Clear everything

% -------------------------------------------------------------- Load File
% Load the ECG file

% Read and ECG signal into memory
load SampleECG.txt;  % Load file into memory
ECGTimeAxis = SampleECG(1:2000,1);   %  copy column 1 which is time axis
ECGVoltage = SampleECG(1:2000,2);   %  copy column 2 which is voltage

fs = 200; % Sampling frequency is fs = 200 Hz
N = length(ECGVoltage); % number of samples
T = 1/fs; % period

% -------------------------------------------------------------- Part 1
% Apply the von Hann lowpass ﬁlter.  Specify the
% ﬁlter in terms of the a and b arrays via the filter command in MATLAB.
% Obtain its freqency responce (magnitude and phase), pole-zero plot, 
% as well as the Fourier spectra of the input and output signals.

% Derive low pass (LP) Filter coefficients
% H(z) = G * (1/4) * [1 + 2z^-1 + z^-2]
% Y(z) = G * (1/4) * X(z) + 2X(x)z^-1 + X(z)z^-2
% y(n) = G * (1/4) * x(n) + 2x(n-1) + x(n-2)
% Apply to form a(1)y(n) = b(1)x(n) + b(2)x(n-1) + b(3)x(n-2)
LPa(1) = 1;
LPb(1) = 1;
LPb(2) = 2;
LPb(3) = 1;

% normalize gain at z=1 to 1
% H(z) = G * [1 + 2z^-1 + z^-2]
% H(z=1) = 1 = G * (1 + 2 + 1)
LPGain = 1/4;
LPb = LPb * LPGain % combine into b coefficients

% Filter the ECG signal
ECGLPFilteredVoltage = filter(LPb, LPa, ECGVoltage);

%Plot the unfiltered ECG Signal vs. Time
LPECGPlot = figure('Name','Von Hann lowpass ﬁlter - ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage);
title('SampleECG1.txt - Unfiltered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis tight;

%Plot the Filtered ECG Signal vs. Time
subplot(2,1,2); plot(ECGTimeAxis, ECGLPFilteredVoltage);
title('SampleECG.txt - Lowpass ﬁltered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis tight;

%Plot the unfiltered ECG Signal vs. Time for 2 Cycles
LPECGPlot2 = figure('Name','Cascade Filter - 2 Cycles ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage);
title('SampleECG.txt - Unfiltered ECG Signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis([0 1.5 -0.5 1.5]);

%Plot the Filtered ECG Signal vs. Time for 2 Cycles
subplot(2,1,2); plot(ECGTimeAxis, ECGLPFilteredVoltage);
title('SampleECG.txt - Cascade filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis([0 1.5 -0.5 1.5]);

% Plot the PSD of the unfiltered signal
LPECGPSDPlot = figure('Name','Von Hann lowpass ﬁlter - PSDs of ECG'); % Create a new figure
subplot(2,1,1); ECG = FindPSD2(ECGVoltage, (length(ECGVoltage)), fs);
title('PSD of unfiltered ECG Signal (Fs = 200 Hz) (N = # Number of samples) (SampleECG.txt)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis([0 (fs/2) -40 5]);

% Plot the PSD of the filtered signal
subplot(2,1,2); ECGLPPSD = FindPSD2(ECGLPFilteredVoltage, (length(ECGLPFilteredVoltage)), fs);
title('Derivative Filter - PSD of Filtered ECG Signal (Fs = 200 Hz) (N = # Number of samples) (SampleECG.txt)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis([0 (fs/2) -40 5]);

% Plot the filter responce
LPFreqResponse = figure('Name','Von Hann lowpass ﬁlter - Frequency & Phase Response'); % Create a new figure
freqz(LPb,LPa, 512, 200);
title('Von Hann lowpass ﬁlter Frequency & Phase Response');

% Plot the filter pole - zero diagram
% H(z) = (Z+1)(Z+1) / Z^2
LPZPlane = figure('Name','Von Hann lowpass ﬁlter - Pole-Zero diagram of the filter'); % Create a new figure
zplane(LPb,LPa)
title('Von Hann lowpass ﬁlter - Pole-Zero diagram of the filter');

% -------------------------------------------------------------- Part 2
% Modify the derivative-based ﬁlter given by Equation 3.47 in the textbook for the
% removal of low-frequency artifacts so that the gain at the maximum frequency present
% in the input signal is unity. Use the filter command and apply the ﬁlter to your
% signal.  Obtain its freqency responce (magnitude and phase), pole-zero plot, 
% as well as the Fourier spectra of the input and output signals.

% Derive derivative filter (DF) coefficients
% H(z) = G * (1/T) * [1 - z^-1] / [1 - 0.995z^-1]
% Assume that T = 1
% [1 - 0.995z^-1]Y(z) = GX(z)[1-z^-1]
% Y(z) - Y(z)0.995z^-1 = GX(z) - GX(z)z^-1
% Y(z) = GX(z) - GX(z)z^-1 + Y(z)0.995Z^-1
% y(n) = Gx(n) - Gx(n-1) + 0.995y(n-1)
% Apply to form a(1)y(n) + a(2)y(n-1) = b(1)x(n) + b(2)x(n-1) + b(3)x(n-2)
DFa(1) = 1;
DFa(2) = -0.995;
DFb(1) = 1;
DFb(2) = -1;

% normalize gain at z=-1 to 1
% H(z) = G (1 - z^-1) / (1 - 0.995z^-1)
% H(z=-1) = 1 = G (1 - (-1) / ( 1 - 0.995(-1)
% 1 = G * (2 / 1.9950)
DFGain = 0.9975;
DFb = DFb * DFGain % combine into b coefficients

% Filter the ECG signal
ECGDerivativeFilteredVoltage = filter(DFb, DFa, ECGVoltage);

%Plot the unfiltered ECG Signal vs. Time
DFECGPlot = figure('Name','Derivative Filter - ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage);
title('SampleECG.txt - Unfiltered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis tight;

%Plot the Filtered ECG Signal vs. Time
% Check out (1:2000 to see some baseline reduction .... need long segment
subplot(2,1,2); plot(ECGTimeAxis, ECGDerivativeFilteredVoltage);
title('SampleECG.txt - Derivative Filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis tight;

%Plot the unfiltered ECG Signal vs. Time for 2 Cycles
DFECGPlot2 = figure('Name','Cascade Filter - 2 Cycles ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage);
title('SampleECG.txt - Unfiltered ECG Signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis([0 1.5 -0.5 1.5]);

%Plot the Filtered ECG Signal vs. Time for 2 Cycles
subplot(2,1,2); plot(ECGTimeAxis, ECGDerivativeFilteredVoltage);
title('SampleECG.txt - Cascade filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis([0 1.5 -0.5 1.5]);

% Plot the PSD of the unfiltered signal
DFECGPSDPlot = figure('Name','Derivative Filter - PSDs of ECG'); % Create a new figure
subplot(2,1,1); ECG = FindPSD2(ECGVoltage, (length(ECGVoltage)), fs);
title('PSD of unfiltered ECG Signal (Fs = 200 Hz) (N = # Number of samples) (SampleECG.txt)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis([0 (fs/2) -40 5]);

% Plot the PSD of the filtered signal
subplot(2,1,2); ECGDFPSD = FindPSD2(ECGDerivativeFilteredVoltage, (length(ECGDerivativeFilteredVoltage)), fs);
title('Derivative Filter - PSD of Filtered ECG Signal (Fs = 200 Hz) (N = # Number of samples) (SampleECG.txt)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis([0 (fs/2) -40 5]);

% Plot the filter responce
DFFreqResponse = figure('Name','Derivative Filter - Frequency & Phase Response'); % Create a new figure
freqz(DFb,DFa, 512, 200);
title('Derivative Filter Frequency & Phase Response');

% Plot the filter pole - zero diagram
DFZPlane = figure('Name','Derivative Filter - Pole-Zero diagram of the filter'); % Create a new figure
zplane(DFb,DFa)
title('Derivative Filter - Pole-Zero diagram of the filter');

% -------------------------------------------------------------- Part 3
% Apply the notch ﬁlter that you designed in Lab 6 for the rejection of 60 Hz to your
% signal. Ensure that the ﬁlter is normalized to have unit gain at DC. Use the filter
% command.  Obtain its freqency responce (magnitude and phase), pole-zero plot, 
% as well as the Fourier spectra of the input and output signals.

% Determine the where the zeros go on the notch filter
OmegaZero = 2 * pi * (60 / 200) % OmegaZero = 2*pi (notch freq / sampling freq)
z1 = cos(OmegaZero) + 1i * sin(OmegaZero)
z2 = cos(-OmegaZero) + 1i * sin(-OmegaZero)

% Derive Notch Filter (NF) coefficients
% Convert the transfer function H(z) to a difference equation y(n)
% H(z) = (z-z1)(z-z2) / z^2 
% H(z) = (Z^2 - (z1+z2)z + z1z2) / z^2
% H(z) = 1 - (z1+z2)z^-1 + (z1z2)z^-2
% H(z) = Y(z) / X(z)
% Y(z) = X(z) * (1 - (z1+z2)z^-1 + (z1z2)z^-2)
% Y(z) = X(z) - X(z)(z1+z2)z^-1 + X(z)(z1z2)z^-2
% y(n) = x(n) - (z1+z2)*x(n-1) + (z1z2)*x(n-2)
% y(n) = b(1)x(n) + b(2)x(n-1) + b(3)x(n-2)
NFa(1) = 1;
NFb(1) = 1;
NFb(2) = -z1 + -z2;
NFb(3) = z1 * z2;

% For low-pass and notch filters we set normzlize the DC gain to 1
% H(z) = Gain * (b(1) + b(2)z^-1 + b(3)z^-2)
% H(z=1) = Gain * (b(1) + b(2)1 + b(3)1) = 1
NFHZEqualOne = NFb(1) + NFb(2) + NFb(3);
NFGain = 1 / NFHZEqualOne
NFb = NFb * NFGain

% Filter the ECG signal
ECGNotchFilteredVoltage = filter(NFb, NFa, ECGVoltage);

%Plot the unfiltered ECG Signal vs. Time
NFECGPlot = figure('Name','Notch Filter - ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage);
title('SampleECG.txt - Unfiltered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis tight;

%Plot the Filtered ECG Signal vs. Time
subplot(2,1,2); plot(ECGTimeAxis, ECGNotchFilteredVoltage);
title('SampleECG.txt - Notch filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis tight;

%Plot the unfiltered ECG Signal vs. Time for 2 Cycles
NFECGPlot2 = figure('Name','Cascade Filter - 2 Cycles ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage);
title('SampleECG.txt - Unfiltered ECG Signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis([0 1.5 -0.5 1.5]);

%Plot the Filtered ECG Signal vs. Time for 2 Cycles
subplot(2,1,2); plot(ECGTimeAxis, ECGNotchFilteredVoltage);
title('SampleECG.txt - Cascade filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis([0 1.5 -0.5 1.5]);

% Plot the PSD of the unfiltered signal
NFECGPSDPlot = figure('Name','Notch Filter - PSDs of ECG'); % Create a new figure
subplot(2,1,1); ECGPSD = FindPSD2(ECGVoltage, (length(ECGVoltage)), fs);
title('PSD of unfiltered ECG Signal (Fs = 200 Hz) (N = # Number of samples) (SampleECG.txt)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis([0 (fs/2) -40 5]);

% Plot the PSD of the filtered signal
subplot(2,1,2); ECGNFPSD = FindPSD2(ECGNotchFilteredVoltage, (length(ECGNotchFilteredVoltage)), fs);
title('Notch Filter - PSD of notch filtered ECG Signal (Fs = 200 Hz) (N = # Number of samples) (SampleECG.txt)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis([0 (fs/2) -40 5]);

% Plot the filter responce
NFFreqResponse = figure('Name','Notch Filter - Frequency & Phase Response'); % Create a new figure
freqz(NFb,NFa, 512, 200);
title('Notch Filter - Frequency & Phase Response');

% Plot the filter pole - zero diagram
NFZPlaneFilter = figure('Name','Notch Filter - Pole-Zero diagram of the filter'); % Create a new figure
zplane(NFb,NFa)
title('Notch Filter - Pole-Zero diagram of the filter');

% -------------------------------------------------------------- Part 4
% Part 4 - Apply all three ﬁlters to the ECG signal in series, and study the combined ﬁlter and
% the result as speciﬁed above.  Obtain its freqency responce (magnitude and phase), 
% pole-zero plot, as well as the Fourier spectra of the input and output signals.

% Find the filter coefficients for the cascade filter
% Hcascade = H(z)lowpass * H(z)derivative & H(z)notchFilter
% Hcascade = LPb * NFb * DFb / LPa * NFa * DFa
Cb = conv(NFb, conv(DFb, LPb))
Ca = conv(NFa, conv(DFa, LPa)) 

% Filter the ECG signal
ECGCascadeFilteredVoltage = filter(Cb, Ca, ECGVoltage);

%Plot the unfiltered ECG Signal vs. Time
CECGPlot = figure('Name','Cascade Filter - ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage);
title('SampleECG.txt - Unfiltered ECG Signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis tight;

%Plot the Filtered ECG Signal vs. Time
subplot(2,1,2); plot(ECGTimeAxis, ECGCascadeFilteredVoltage);
title('SampleECG.txt - Cascade filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis tight;

%Plot the unfiltered ECG Signal vs. Time for 2 Cycles
CECGPlot2 = figure('Name','Cascade Filter - 2 Cycles ECG signal'); % Create a new figure
subplot(2,1,1); plot(ECGTimeAxis, ECGVoltage);
title('SampleECG.txt - Unfiltered ECG Signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis([0 1.5 -0.5 1.5]);

%Plot the Filtered ECG Signal vs. Time for 2 Cycles
subplot(2,1,2); plot(ECGTimeAxis, ECGCascadeFilteredVoltage);
title('SampleECG.txt - Cascade filtered ECG signal');
xlabel('Time (Sec)'); ylabel('ECG (mV)'); axis([0 1.5 -0.5 1.5]);

% Plot the PSD of the unfiltered signal
CECGPSDPlot = figure('Name','Cascade Filter - PSDs of ECG'); % Create a new figure
subplot(2,1,1); ECGPSD = FindPSD2(ECGVoltage, (length(ECGVoltage)), fs);
title('PSD of unfiltered ECG Signal (Fs = 200 Hz) (N = # Number of samples) (SampleECG.txt)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis([0 (fs/2) -40 5]);

% Plot the PSD of the filtered signal
subplot(2,1,2); ECGCPSD = FindPSD2(ECGCascadeFilteredVoltage, (length(ECGCascadeFilteredVoltage)), fs);
title('PSD of cascade filtered ECG Signal (Fs = 200 Hz) (N = # Number of samples) (SampleECG.txt)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis([0 (fs/2) -40 5]);

% Plot the filter responce
CFreqResponse = figure('Name','Cascade Filter - Frequency & Phase Response'); % Create a new figure
freqz(Cb,Ca, 512, 200);
title('Cascade Filter - Frequency & Phase Response');

% Plot the filter pole - zero diagram
CZPlaneFilter = figure('Name','Cascade Filter - Pole-Zero diagram of the filter'); % Create a new figure
zplane(Cb,Ca)
title('Cascade Filter - Pole-Zero diagram of the filter');
