function PSD = FindPSD(SignalOfInterest, N, fs)
%FindPSD - Returns a PSD of a signal.
% We will pass the signal of interest, number of samples (N) and the
% sampling frequency (fs) to this function.  The function returns a
% structure containing the PSD of the signal.
% Inputs:
% SignalOfInterest = Signal or segment of interest
% N = Number of points
% fs = sample rate of Signal

% ---------------------------------------------------- Obtain the PSD
% Select the segment of interest

% Subtract the mean of the segment; this will set the DC component in the PSD to zero.
% Otherwise the DC component will be above -100 dB.
SignalOfInterest = SignalOfInterest - mean(SignalOfInterest);

% Obtain the Fourier transform (FT) of the segment. The Matlab command for the com-
%putation of the FT is fft (representing the fast Fourier transform or FFT algorithm).
% Obtain the squared magnitude of the FT.
PSD = ( abs( fft(SignalOfInterest, N) ) ).^2; % magnitude squared

% Select the first half of the PSD above. Prepare the corresponding frequency axis in
% Hz, spanning the range [0, fs/2], where fs is the sampling frequency of the signal.
% PSD is symmetric around fs/2 so we only need to show half of it.

PSDAxis = (0 : (ceil(length(PSD)) - 1 )) * (fs / length(PSD) );
%Plot and analyze the result. Ensure that the axes are labeled in Hz and dB.
%Compare signals & bandwidth

plot(PSDAxis, PSD);
title('Sample Segment PSD');
xlabel('Freq (Hz)'); ylabel('Magnitude'); axis auto;

end