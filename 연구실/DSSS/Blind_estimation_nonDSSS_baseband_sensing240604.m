clc; 
clear;
close all

N_bit=1000;
%임의의 PSK 신호 생성
SNR=-10;
M=2;
bits=2*randi([0,1], N_bit,1)-1; %BPSK
Tb=1/(10^5); %bit duration
oversamplingrate=4; %1chip=4samples
fc=1*10^7; %80MHz=8*10^7Hz
f1=fc*0.95;
f2=fc*1.2;
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
    'OutputSamplesPerSymbol', oversamplingrate*length(pn), ...
    'FilterSpanInSymbols', 10); %OutputSamplesPerSymbol 사실상 chip당 샘플수
waveform0=rcfilter(bits).';
Tx=waveform0.*cos(2*pi*fc*t+pi/6);
Rx=awgn(waveform0, SNR, 'measured');%SNR db scale

f=linspace(-fs/2,fs/2-fs/length(Tx),length(Tx));

%신호 plot
figure(1)
subplot(3,1,1)
plot(t, oversampled_bits)
title("Original signal")
xlim([0, 1/(10^5)])
subplot(3,1,2)
plot(t, oversampled_spreaded_bits)
title("Spreaded signal")
xlim([0, 1/(10^5)])
subplot(3,1,3)
plot(t, waveform0)
title("Spreaded signal in baseband")
xlim([0, 1/(10^6)])
set(gcf, 'Color', 'w'); % figure 배경을 흰색으로 설정
set(gca, 'Color', 'w'); % axes 배경을 흰색으로 설정



ts=1; %cyclic autocorrelation tau
figure(2)
[cf, tau, Tx_psd, p] = cyclic_autocorr(Tx, ts, fs);
%index4=find(cf == 0); %cf=0의 spectral line삭제
%Rx_psd(index4) = 0;
plot(cf, Tx_psd);
xlabel('x=cylic frequency')
ylabel('y=cyclic autocorrelation')
set(gcf, 'Color', 'w'); % figure 배경을 흰색으로 설정
set(gca, 'Color', 'w'); % axes 배경을 흰색으로 설정

[sortedValues, sortedIndices] = sort(Tx_psd, 'descend');
top3=sortedIndices(1:3); 
for k=1:3
    if cf(sortedIndices(k)) >=1
        estimated_fc=cf(sortedIndices(k))/2;
    end
end
sprintf("DSSS signal frequency: %d", fc)
sprintf("Estimated DSSS signal frequency: %d", estimated_fc)

%중심주파수 fc를 통해서 down conversion 되어있다고 가정하고 baseband에서 CAF

figure(3)
waveform0_FFT=abs(fftshift(fft(waveform0)));
plot(f, waveform0_FFT)
title("baseband waveform")


figure(4)
[cf, tau, Rx_psd, p] = cyclic_autocorr(Rx, ts, fs);
plot(cf, Rx_psd);
xlabel('x=cylic frequency')
ylabel('y=cyclic autocorrelation')
set(gcf, 'Color', 'w'); % figure 배경을 흰색으로 설정
set(gca, 'Color', 'w'); % axes 배경을 흰색으로 설정


Rx_psd0=Rx_psd;
%CFAR 
index0=find(cf == 0); %cf=0의 spectral line삭제
Rx_psd0(index0) = 0;
cells=length(Rx_psd0)/(10^4)/2;
numGuardCells=round(cells*0.01);
numRefCells=round(cells*0.99);
thresholdFactor=5;
[cfar_targets] = cfar_ca_1D(Rx_psd0, numGuardCells, numRefCells, thresholdFactor);
figure(5)
plot(cf, Rx_psd0);
xlabel('x=cylic frequency')
ylabel('y=cyclic autocorrelation')
set(gcf, 'Color', 'w'); % figure 배경을 흰색으로 설정
set(gca, 'Color', 'w'); % axes 배경을 흰색으로 설정


%피크 탐지
hold on;
m=mean(Rx_psd0);
peaks=(cfar_targets.') .* Rx_psd0;
index1=find((peaks > m));
[chiprate, indexi] = max(cf(index1));
scatter(cf(index1), Rx_psd0(index1), 'r*');
title('CFAR Target Detection');
legend('Input Signal', 'Detected Peaks');
sprintf("Estimated chip rate: %d", chiprate)

