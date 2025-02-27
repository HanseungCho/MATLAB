clc; 
clear;
close all
N_bit=3000;
%임의의 PSK 신호 생성

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
fs=(1/Tb)*length(pn)*oversamplingrate;%sampling frequency = 10^9Hz

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
Tx=waveform0.*cos(2*pi*fc*t+0);
Rx=awgn(Tx, 0, 'measured');
downRx=Tx.*cos(2*pi*fc*t);

f=linspace(-fs/2,fs/2,length(Tx));
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
waveform0_FFT=abs(fftshift(fft(waveform0)));
Tx_FFT=abs(fftshift(fft(Tx)));
Rx_FFT=abs(fftshift(fft(Rx)));
downRx_FFT=abs(fftshift(fft(downRx)));

figure(3)
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
title("passband waveform")
ylim([0 3*10^4])
subplot(1,4,4)
plot(f, downRx_FFT, 'b')
title("baseband waveform downconversed")
ylim([0 3*10^4])

sprintf("DSSS signal frequency: %d, signal1 frequency: %d, signal2_frequency: %d", fc, f1, f2)

ts=1;
figure(6)
%Cylic autocorrelation
[cf, tau, Tx_psd, p] = cyclic_autocorr(waveform0, ts, fs);
%index0=find(cf == 0);%cf=0의 spectral line삭제
%Tx_psd(index0) = 0;
plot(cf, Tx_psd);
xlabel('x=cylic frequency')
ylabel(sprintf('y=DSSS only cyclic autocorrelation (tau = %d)', tau(p)*(1/fs)))

figure(8)
[cf, tau, Rx_psd, p] = cyclic_autocorr(Rx, ts, fs);
%index4=find(cf == 0);%cf=0의 spectral line삭제
%Rx_psd(index4) = 0;
plot(cf, Rx_psd);
xlabel('x=cylic frequency')
ylabel(sprintf('y= Tx,Tx_s1,2,noise cyclic autocorrelation (tau = %d)', tau(p)*(1/fs)))

[tau, p, quadratic] = quadratic(downRx, ts);
%figure(9)
%[coeffs, a, b] = customHaarCWT(quadratic, Tb, fs);
%imagesc(b, a/fs, abs(coeffs));
%xlabel('Time'); ylabel('Scale');
%title('Haar Wavelet Transform Coefficients');
%colorbar;
%[maxvalue, linearindex]=max(abs(coeffs(:)));
%[rowindex, columnindex]=ind2sub(size(abs(coeffs(:))), linearindex);

%figure(9)
%[cfs, frq] =cwt(quadratic);
%surface(t, frq, abs(cfs));
