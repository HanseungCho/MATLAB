clc; clear;
%%
cfgSU = wlanHESUConfig;
cfgSU.ExtendedRange = false;      % Do not use extended-range format
cfgSU.ChannelBandwidth = 'CBW20'; % Channel bandwidth
cfgSU.APEPLength = 1000;          % Payload length in bytes
cfgSU.MCS = 0;                    % Modulation and coding scheme
cfgSU.ChannelCoding = 'LDPC';     % Channel coding
cfgSU.NumSpaceTimeStreams = 1;    % Number of space-time streams
cfgSU.NumTransmitAntennas = 1;    % Number of transmit antennas
%%
psdu = randi([0 1],getPSDULength(cfgSU)*8,1,'int8'); % Random PSDU
txSUWaveform = wlanWaveformGenerator(psdu,cfgSU);    % Create packet
%%
cfgExtSU = cfgSU;
cfgExtSU.ExtendedRange = true;  % Enable extended-range format
cfgExtSU.Upper106ToneRU = true; % Use only upper 106-tone RU

% Generate a packet
psdu = randi([0 1],getPSDULength(cfgExtSU)*8,1,'int8'); % Random PSDU
txExtSUWaveform = wlanWaveformGenerator(psdu,cfgExtSU);   % Create packet
%%
fs = wlanSampleRate(cfgExtSU); % Get baseband sample rate
ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgExtSU);
fftsize = ofdmInfo.FFTLength; % Use the data field fft size
rbw = fs/fftsize; % Resoluton bandwidth
spectrumScope = spectrumAnalyzer(SampleRate=fs,...
    RBWSource='property',RBW=rbw,...
    AveragingMethod='exponential',ForgettingFactor=0.25,...
    YLimits=[-50,20],...
    Title='HE Extended-Range SU with Active Upper 106-Tone RU');
spectrumScope.ViewType = 'Spectrum and Spectrogram';
spectrumScope.TimeSpanSource = 'Property';
spectrumScope.TimeSpan = length(txExtSUWaveform)/fs;
spectrumScope(txExtSUWaveform)
%%
figure;
ind = wlanFieldIndices(cfgExtSU);
t = (0:(double(ind.LLTF(2))-1))/fs*1e6;
plot(t,20*log10(movmean(abs(txSUWaveform(1:ind.LLTF(2))),20)),'-b')
hold on;
plot(t,20*log10(movmean(abs(txExtSUWaveform(1:ind.LLTF(2))),20)),'-r')
grid on;
title('Power of L-STF and L-LTF (1 us Moving Average)');
xlabel('Time (us)');
ylabel('Power (dBW)');
legend('HE SU','HE Extended-Range SU','Location','SouthWest');
%%
