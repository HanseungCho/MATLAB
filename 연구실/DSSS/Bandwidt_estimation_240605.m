clc; clear; close all;

% Parameters
M = 2;
in = randi([0, 1], 60000, 1);
usedInput = in(1:log2(M) * (floor(length(in) / log2(M))));
symbolInput = bit2int(usedInput, log2(M));
T = 1 / 1000;
Rs = 1 / T;
tx = pskmod(symbolInput, M, 0);
oversamplingRate = 4;
fs = oversamplingRate / T;

% Raised Cosine Filter
rcFilter = comm.RaisedCosineTransmitFilter('Shape', 'Square root', ...
    'RolloffFactor', 0.2, ...
    'OutputSamplesPerSymbol', oversamplingRate, ...
    'FilterSpanInSymbols', 10);
waveform = rcFilter(tx);

% Add noise
rx = awgn(waveform, 10, 'measured');

% Estimate Welch Spectrum
[pxx, f] = pwelch(rx, 1024, [], [], fs);

% Plot Welch Spectrum
figure;
f=f-fs/2;
plot(f, 10*log10(fftshift(pxx)));
title('Welch Power Spectral Density Estimate');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

% Cumulative Sum of Spectrum
cumSumPxx = cumsum(fftshift(pxx));
figure;
plot(f, cumSumPxx)
title('Cumulative sum sum of Welch Power Spectral Density Estimate');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

%% 가우시안 커널 생성
%sigma = 2;  % 표준 편차
%kernel_size = 6 * sigma;  % 커널의 크기 (표준 편차의 6배 범위 내에서 자름)
%x = -kernel_size:kernel_size;
%gaussian_kernel = (1/(sqrt(2*pi)*sigma)) * exp(-x.^2/(2*sigma^2));
%
%% 커널 정규화
%gaussian_kernel = gaussian_kernel / sum(gaussian_kernel);
%
%cumSumPxx = conv(cumSumPxx, gaussian_kernel, 'same');

% Second derivative
second_derivative = diff(diff(cumSumPxx));
% Find inflection points (zero-crossings of the second derivative)
zero_crossings = find(second_derivative(1:end-1) .* second_derivative(2:end) < 0);

% Estimate bandwidth as the difference between first and last zero-crossing
bandwidth_idx = zero_crossings(end) - zero_crossings(1);
bandwidth = f(bandwidth_idx);

disp(['Estimated Bandwidth: ', num2str(bandwidth), ' Hz']);

% Plot the inflection points on the cumulative sum plot
figure;
plot(f, cumSumPxx, '-b'); hold on;
plot(f(zero_crossings), cumSumPxx(zero_crossings), 'ro');
title('Cumulative sum of Welch Power Spectral Density Estimate with Inflection Points');
xlabel('Frequency (Hz)');
ylabel('Cumulative Power');
grid on;
