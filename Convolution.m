% 2.1 Convolution
% Convolves two signals x(n) and h(n) and plots original signals and result
% of convolution

clear all % Clear all Varables
close all % Clear all windows

X = [ 0 2 3 4 -5]; % x(n)
H = [0.5 0.5]; % h(n)
Y = conv(X, H); % y(n) = x(n) Convolved with h(n)

figure('Name','Convolution of x(n) with h(n)')

subplot(3,1,1), stem(X) % Plot x(n)
title('x(n)');
xlabel('n');
ylabel('Amplitude');
axis([0 7 -6 6]); % Setup axis so it lines up with y(n)'s axis

subplot(3,1,2), stem(H) % Plot h(n)
title('h(n)')
xlabel('n');
ylabel('Amplitude');
axis([0 7 0 1]); % Setup axis so it lines up with y(n)'s axis

subplot(3,1,3), stem(Y)
title('y(n) = Convolution of x(n) with h(n)') % Plot y(n)
xlabel('n');
ylabel('Amplitude');
axis([0 7 -6 6]);