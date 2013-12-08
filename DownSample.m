function DownSampledECG = DownSample(ECG, Fs)
% FindQRS - Low pass filters the signal, down samples it to 200 Hz
% Desc

% Preﬁlter the ECG channel only using a second-order Butterworth lowpass 
% ﬁlter with a cutoﬀ frequency of 60 Hz
[b,a] = butter(2, (60 / (Fs/2)));
LowPassECG = filter(b, a, ECG);

% downsample by a factor of ﬁve before applying the Pan–Tompkins algorithm. 
DownSampledECG = LowPassECG(5:5:end); % takes every fifth sample