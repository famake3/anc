pkg load signal
% test of math principles

Fs = 44100;
arrival_a = 0.0*Fs
arrival_b = 0.5*Fs

x = -Fs*2:Fs*2;
t = x / Fs;
N = size(x)(2)
a_signal =     normpdf(x, arrival_a, Fs*0.2) * Fs*0.5;
a_background = normpdf(x, -Fs*1.5, Fs*0.2) * Fs*0.5;
a = a_signal + a_background;

b_signal =     normpdf(x, arrival_b, Fs*0.2) * Fs*0.5;
%b_background = normpdf(x, -Fs*1.5, Fs*0.2) * Fs*0.5;
b = b_signal ; %+ b_background;

figure;
title("Original signals");
subplot(2, 1, 1);
plot(t, a)
subplot(2, 1, 2);
plot(t, b)

% Length of the FFT:
% Should have (2*N - 1) for 
ft_a = fft(a, 2*N - 1);
ft_b = fft(b, 2*N - 1);

ftprod = ft_a .* conj(ft_b);
corr_scrambled = ifft(ftprod);
corr = [corr_scrambled(N+1:end) corr_scrambled(1:N)];
lag = ( -(N-1):(N-1) ) / Fs;
% Meaning of the correlation function and the lag:
% Data at lag l represents how much the first signal (A) 
% lags behind the second. How much greater the arrival time
% is in A vs. B.

figure;
title("Correlation function");
plot(lag, corr)


% Try to only keep data which are correlated

figure;
title("Correlation function");
plot(lag, corr)
