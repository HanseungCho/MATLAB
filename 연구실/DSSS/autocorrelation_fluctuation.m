clc; 
clear;
close all

N_bit=10000;
%임의의 PSK 신호 생성
SNR=-5;

bits=2*randi([0,1], N_bit,1)-1; %BPSK
Tb=1/(10^5); %bit duration
oversamplingrate=4; %1chip=4samples
fc=1*10^7; %80MHz=8*10^7Hz
PN=comm.PNSequence('Polynomial',[1 0 0 0 1 1 0 1], 'SamplesPerFrame', 127, 'InitialConditions',[0 0 0 0 0 0 1]);
pn=PN();
Processing_Gain=length(pn);
Rc=Processing_Gain/Tb;%chip rate => chip rate
fs=(1/Tb)*length(pn)*oversamplingrate;

for k=1:length(bits)
    over_bits(length(pn)*(k-1)+1:length(pn)*k)=bits(k);
    spreaded_bits(length(pn)*(k-1)+1:length(pn)*k)=bits(k)*(2*pn-1);
end
%펄스로 나타낸 심볼, 칩
oversampled_bits=repelem(over_bits,oversamplingrate);
oversampled_spreaded_bits=repelem(spreaded_bits,oversamplingrate);
phase=rand(1)*2*pi;
t=linspace(0,(length(bits)*Tb)-(1/fs),length(oversampled_spreaded_bits));

%실제로는 rrc필터를 형태로 신호를 생성하게됨
rcfilter = comm.RaisedCosineTransmitFilter('Shape', 'Square root', ...
    'RolloffFactor', 0.2, ...
    'OutputSamplesPerSymbol', oversamplingrate, ...
    'FilterSpanInSymbols', 10); %OutputSamplesPerSymbol 사실상 chip당 샘플수
waveform0=rcfilter(spreaded_bits.').';
Tx=waveform0;
Rx=awgn(Tx, SNR, 'measured');%SNR db scale
Rx=circshift(Rx,round((length(Rx)/N_bit)*rand(1))); %DSSS 신호에 랜덤한 timing offset 부여

%윈도우 생성
M=(length(Rx)/N_bit)*10;
window=zeros(length(Rx)/M, M);
T=length(Rx)/M;
for i=1:M
    window(i,1:T)=Rx(T*(i-1)+1: T*i);
    [c, lags]=xcorr(window(i,:));
    ACM(i,:)=c; %Autocorrelation matrix
end
g=mean(abs(ACM).^2);
plot(lags*(1/fs),g)
title("DSSS signal")

