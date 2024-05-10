clc; 
clear;
close all
for p
N_bit=5000;
%임의의 PSK 신호 생성
SNR=0;
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
    'OutputSamplesPerSymbol', oversamplingrate, ...
    'FilterSpanInSymbols', 10); %OutputSamplesPerSymbol 사실상 chip당 샘플수
waveform0=rcfilter(spreaded_bits.').';
Tx=waveform0.*cos(2*pi*fc*t+pi/6);
Rx=awgn(Tx, SNR, 'measured');%SNR db scale


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



ts=1; %cyclic autocorrelation tau
figure(2)
[cf, tau, Rx_psd, p] = cyclic_autocorr(Rx, ts, fs);
%index4=find(cf == 0); %cf=0의 spectral line삭제
%Rx_psd(index4) = 0;
plot(cf, Rx_psd);
xlabel('x=cylic frequency')
ylabel('y=cyclic autocorrelation')
[sortedValues, sortedIndices] = sort(Rx_psd, 'descend');
top3=sortedIndices(1:3); 
for k=1:3
    if cf(sortedIndices(k)) >=1
        estimated_fc=cf(sortedIndices(k))/2;
    end
end
sprintf("DSSS signal frequency: %d", fc)
sprintf("Estimated DSSS signal frequency: %d", estimated_fc)

%down conversion
cosmix_Rx=Rx.*cos(2*pi*estimated_fc*t);
[vl,indexh]=min(abs(f-estimated_fc));
[vh,indexl]=min(abs(f-(-estimated_fc)));
downRx=cosmix_Rx;
downRx_FFT=fftshift(fft(downRx));
downRx_FFT(1:indexl)=0;
downRx_FFT(indexh:end)=0;


figure(3)
waveform0_FFT=abs(fftshift(fft(waveform0)));
Tx_FFT=abs(fftshift(fft(Tx)));
Rx_FFT=abs(fftshift(fft(Rx)));
subplot(1,4,1)
plot(f, waveform0_FFT)
title("baseband waveform")
ylim([0 3*10^4])
subplot(1,4,2)
plot(f, Tx_FFT)
title("passband waveform")
ylim([0 3*10^4])
subplot(1,4,3)
plot(f, Rx_FFT)
title("passband waveform with noise")
ylim([0 3*10^4])
subplot(1,4,4)
plot(f, abs(downRx_FFT), 'b')
title("baseband waveform downconversed")
ylim([0 3*10^4])

figure(4)
basebandRx=ifft(ifftshift(downRx_FFT));
[cf, tau, basebandRx_psd, p] = cyclic_autocorr(basebandRx, ts, fs);
plot(cf, basebandRx_psd);
xlabel('x=cylic frequency')
ylabel('y=cyclic autocorrelation')

basebandRx_psd0=basebandRx_psd;
%CFAR 
index0=find(cf == 0); %cf=0의 spectral line삭제
basebandRx_psd0(index0) = 0;
cells=length(basebandRx_psd0)/(10^4)/2;
numGuardCells=round(cells*0.01);
numRefCells=round(cells*0.99);
thresholdFactor=10;
[cfar_targets] = cfar_ca_1D(basebandRx_psd0, numGuardCells, numRefCells, thresholdFactor);
figure(5)
plot(cf, basebandRx_psd0, 'b');
hold on;
m=mean(basebandRx_psd0);
peaks=(cfar_targets.') .* basebandRx_psd0;
index1=find((peaks > m));
[chiprate, indexi] = max(cf(index1));
scatter(cf(index1), basebandRx_psd0(index1), 'r*');
title('CFAR Target Detection');
legend('Input Signal', 'Detected Peaks');
sprintf("Estimated chip rate: %d", chiprate)
%평균 보다 작아지는 지점 찾기 => 해당 지점을 평균을로 하는 gaussian kernel을 곱함
[index1v, index1index]=max(basebandRx_psd0<m);
