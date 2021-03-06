function PSD = FindPSD2(SignalOfInterest, N, fs)
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
FT_Mag = ( abs(fft(SignalOfInterest, N)) ).^2; % magnitude squared

% Set the frequency axis (scale x in Hz)
FTAxis = (0:length(FT_Mag) - 1) * ( fs / length(FT_Mag));

% Normalize the squared magnitude spectrum by dividing by its maximum.
% Convert the result to dB: take the logarithm to base ten of the result above and multiply
% by ten. This is an estimate of the PSD of the original signal in dB.
% Only multiplied by 10 because we squared it, if not use 20

PSD = 10 * log10 (FT_Mag / max(FT_Mag));

% Select the first half of the PSD above. Prepare the corresponding frequency axis in
% Hz, spanning the range [0, fs/2], where fs is the sampling frequency of the signal.
% PSD is symmetric around fs/2 so we only need to show half of it.

PSDAxis = (0 : (ceil(length(FT_Mag) / 2) - 1 )) * (fs / length(FT_Mag) );
PSD = PSD(1:length(PSDAxis));

%Plot and analyze the result. Ensure that the axes are labeled in Hz and dB.
%Compare signals & bandwidth

plot(PSDAxis, PSD);

end