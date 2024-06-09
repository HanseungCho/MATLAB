clc;clear;close all
M=4;
in=randi([0 1], 60000, 1);
usedInput = in(1:log2(M)*(floor(length(in)/log2(M))));
symbolInput = bit2int(usedInput, log2(M));
T=1/(1000);
Rs=1/T;
waveform=pskmod(symbolInput, M, 0);
oversamplingrate=4;
fs=oversamplingrate/T;

rcfilter = comm.RaisedCosineTransmitFilter('Shape', 'Square root', ...
            'RolloffFactor', 0.22, ...
            'OutputSamplesPerSymbol', oversamplingrate, ...
            'FilterSpanInSymbols', 10); %OutputSamplesPerSymbol 사실상 chip당 샘플수
tx=rcfilter(waveform);
rx=awgn(tx, 0, 'measured');
psd=zeros(length(rx), 1);
for tau=0:100
    [f, mem] = cyclic_autocorr(rx, tau, fs);
    psd=psd+mem;
end
psd=psd/(tau+1);
plot(f, db(1000*psd/length(psd)))
ylabel("dBm")
xlim([f(find(f ==0)+10) f(end)])
