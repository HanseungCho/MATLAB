clc; 
clear;
close all

N_bit=100;
%임의의 PSK 신호 생성
SNR=-100;
M=2;
bits=2*randi([0,1], N_bit,1)-1; %BPSK
Tb=1/(10^5); %bit duration
oversamplingrate=8; %1chip=8samples
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
%Rx=Tx;

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
[cf, tau, Rx_psd, p] = cyclic_autocorr(Rx, ts, fs);
%index4=find(cf == 0); %cf=0의 spectral line삭제
%Rx_psd(index4) = 0;
plot(cf, Rx_psd);
xlabel('x=cylic frequency')
ylabel('y=cyclic autocorrelation')
set(gcf, 'Color', 'w'); % figure 배경을 흰색으로 설정
set(gca, 'Color', 'w'); % axes 배경을 흰색으로 설정

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
%title("signal with noise")
xlabel('x=frequency')
ylabel('y=FFT')
ylim([0 3*10^4])
subplot(1,4,4)
plot(f, abs(downRx_FFT), 'b')
title("baseband waveform downconversed")
ylim([0 3*10^4])
set(gcf, 'Color', 'w'); % figure 배경을 흰색으로 설정
set(gca, 'Color', 'w'); % axes 배경을 흰색으로 설정


figure(4)
basebandRx=ifft(ifftshift(downRx_FFT));
[cf, tau, basebandRx_psd, p] = cyclic_autocorr(basebandRx, ts, fs);
plot(cf, basebandRx_psd);
xlabel('x=cylic frequency')
ylabel('y=cyclic autocorrelation')
set(gcf, 'Color', 'w'); % figure 배경을 흰색으로 설정
set(gca, 'Color', 'w'); % axes 배경을 흰색으로 설정
%%%%%%%%%%%%여기까지가 기저대역 신호 생성%%%%%%%%%%%%
CFAR=1/(10^5);
sample_shift=0:1; %한 심볼 내에서 tau를 참고
L=length(basebandRx)/4+1;
w = kaiser(L,1);
spectral_window=fft(w);
detection=zeros(1,length(cf)-L);
tic;
for k=1:length(sample_shift)
    shifted_input=circshift(basebandRx, sample_shift(k)); %sample_shift는 벡터
    quadratic=basebandRx.*shifted_input; %x(t)x*(t-tau)
    Fm(k,:)=fftshift(fft(quadratic));
    cvm(k,:)=(1/length(quadratic))*Fm(k,:);  %covarinace vector=[real(cvm(:,i)) imag(cvm(:,i))]   
end
o=find(cf == 15360000);
dk=find(cf == 3700000);
%alpha=1+(L-1)/2:length(cf)-((L-1)/2)alpha=1+(L-1)/2:length(cf)-((L-1)/2)
for alpha=o:o
    S = zeros(length(sample_shift), length(sample_shift));
    Scj = zeros(length(sample_shift), length(sample_shift));
    for m = 1 : length(sample_shift)
        for n = 1 : length(sample_shift)
            fn=Fm(n,:);
            fm=Fm(m,:);
            fv1=flip(fn(alpha-(L-1)/2:alpha+(L-1)/2));
            fv2=fm(alpha-(L-1)/2:alpha+(L-1)/2);
            fv3=conj(fn(alpha-(L-1)/2:alpha+(L-1)/2));
            
            S(m, n)=fv1.*fv2*spectral_window; % 윈도우 함수의 인덱스가 음수가 되지 않도록 조정
            Scj(m, n)=fv2.*fv3*spectral_window;
        end
    end
    S = S / (length(basebandRx) * L);
    Scj = Scj / (length(basebandRx) * L);
    sig1=real((S+Scj)/2);
    sig2=imag((S-Scj)/2);
    sig3=imag((S+Scj)/2);
    sig4=real((Scj-S)/2);
    CM=[sig1 sig2;sig3 sig4];
    v=[real(cvm(:,alpha).') imag(cvm(:,alpha).')];
    ML=abs(length(basebandRx)*v*inv(CM)*(v.'));
    syms gam
    threshold=chi2inv(1-CFAR,length(v)); %자유도: Covariance vector 길이
    if ML >= threshold
        detection(alpha-(L-1)/2)=1;
    end
end
if sum(detection) > 0
    disp("DSSS signal detection!")
else
    disp("No siganl")
end
toc;
