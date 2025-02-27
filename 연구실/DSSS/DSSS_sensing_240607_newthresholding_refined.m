clc; 
clear;
close all
%cylic covariance(fft)한 것을 fftshift하지 않고 periodic function이라고 생각하고 
N_bit=80;
%임의의 PSK 신호 생성
SNR=-10:0;
M=2;

Tb=1/(10^5); %bit duration
oversamplingrate=8; %1chip=8samples
PN=comm.PNSequence('Polynomial',[1 1 1 1 0 1], 'SamplesPerFrame', 31, 'InitialConditions',[0 0 0 0 1]);
pn=PN();
Processing_Gain=length(pn);
Rc=Processing_Gain/Tb;%chip rate => chip rate
fs=(1/Tb)*length(pn)*oversamplingrate;
CFAR=1/(10^5);
sample_shift=0:0; %한 심볼 내에서 tau를 참고

simulation=10^3;
time=zeros(length(SNR), simulation);
DP=zeros(1,length(SNR));
for u=1:length(SNR)
    sensing=0;
    for q=1:simulation
        bits=2*randi([0,1], N_bit,1)-1; %BPSK
        pn=circshift(pn, randi([0 length(pn)-1],1));
        for k=1:length(bits)
            over_bits(length(pn)*(k-1)+1:length(pn)*k)=bits(k);
            spreaded_bits(length(pn)*(k-1)+1:length(pn)*k)=bits(k)*(2*pn-1);
        end
          
        %실제로는 rrc필터를 형태로 신호를 생성하게됨
        rcfilter = comm.RaisedCosineTransmitFilter('Shape', 'Square root', ...
            'RolloffFactor', 0.22, ...
            'OutputSamplesPerSymbol', oversamplingrate, ...
            'FilterSpanInSymbols', 10); %OutputSamplesPerSymbol 사실상 chip당 샘플수
        
        waveform0=rcfilter(spreaded_bits.').';
        Rx0=awgn(waveform0, SNR(u), 'measured');%SNR db scale
        Rx=Rx0(1:2^14); %FFT 연산의 속도를 고려하여 2의 power로
        %%%%%%%%%%%%여기까지가 기저대역 신호 생성%%%%%%%%%%%%
        cf=linspace(-fs/2, fs/2-fs/length(Rx), length(Rx));    
            
        tic;
        for k=1:length(sample_shift)
            shifted_input=circshift(Rx, sample_shift(k)); %sample_shift는 벡터
            quadratic=Rx.*shifted_input; %x(t)x*(t-tau)
            Fm(k,:)=fftshift(fft(quadratic)); 
        end
        cvm=Fm/length(quadratic);
        o=find(cf == 0);
        L=length(Rx)/4+1;
        w = kaiser(L,1);
        cvm(:,o)=0;
        [minc, maxc]=bounds(abs(cvm(1,o+1:length(cvm))));
        indd=find(abs(cvm) == maxc);
        %plot(abs(cvm))
        lo=(minc+maxc)/2;
        ind0=find(cvm(1,:) > lo);
        ind=ind0((ind0>o)&(ind0<(length(cf)-((L-1)/2))));
        Xa=abs(cvm(1,ind))+cf(ind).*abs(cvm(1,ind));
        [Xav, Xai]= sort(Xa, 'descend');
        if length(Xai) >= 10
            top10Xa=Xai(1:10);
            top10index=ind(top10Xa);
        else
            top10Xa=Xai(1:length(Xai));
            top10index=ind(top10Xa);
        end
        detection=zeros(1,length(cf)-L);
        %alpha=1+(L-1)/2:length(cf)-((L-1)/2) two-sided range => 그런데 양수쪽만
        %봐도 상관 없다. 
        
        for alpha=1:length(top10index)
            S = zeros(length(sample_shift), length(sample_shift));
            Scj = zeros(length(sample_shift), length(sample_shift));
            fv1=flip(Fm(top10index(alpha)-(L-1)/2:top10index(alpha)+(L-1)/2));
            fv2=Fm(top10index(alpha)-(L-1)/2:top10index(alpha)+(L-1)/2);
            fv3=conj(Fm(top10index(alpha)-(L-1)/2:top10index(alpha)+(L-1)/2));
            
    
            S=fv1.*fv2*w; % 윈도우 함수의 인덱스가 음수가 되지 않도록 조정
            Scj=fv2.*fv3*w;
            
            S = S / (length(Rx) * L);
            Scj = Scj / (length(Rx) * L);
            sig1=real((S+Scj)/2);
            sig2=imag((S-Scj)/2);
            sig3=imag((S+Scj)/2);
            sig4=real((Scj-S)/2);
            CM=[sig1 sig2;sig3 sig4];
            v=[real(cvm(:,top10index(alpha)).') imag(cvm(:,top10index(alpha)).')];
            ML=length(Rx)*v/CM*(v.');
            
            threshold=chi2inv(1-CFAR,length(v)); %자유도: Covariance vector 길이
            if ML >= threshold
                detection(alpha)=1;
            end
        end
        sensing=sensing+sign(sum(detection));
        time(u,q)=toc;
    end
    DP(u)=sensing/simulation;
end
figure(1)
plot(SNR, DP)
title("PN=31, CFAR=1/(10^5), tau=0:0, cyclic covariance mean, Number of sample=2^14")
xlabel('SNR(dB)')
ylabel('Detection Probability')