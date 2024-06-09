clc; 
clear;
close all

N_bit=10;
%임의의 PSK 신호 생성
SNR=1:1;
M=2;

Tb=1/(10^5); %bit duration
oversamplingrate=8; %1chip=8samples
PN=comm.PNSequence('Polynomial',[1 0 0 0 1 1 0 1], 'SamplesPerFrame', 127, 'InitialConditions',[0 0 0 0 0 0 1]);
pn=PN();
Processing_Gain=length(pn);
Rc=Processing_Gain/Tb;%chip rate => chip rate
fs=(1/Tb)*length(pn)*oversamplingrate;
CFAR=1/(10^5);
sample_shift=0:0; %한 심볼 내에서 tau를 참고


simulation=100;
DP=zeros(1,length(SNR));
for u=1:length(SNR)
    sensing=0;
    for q=1:simulation
        bits=2*randi([0,1], N_bit,1)-1; %BPSK
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
        Rx=awgn(waveform0, SNR(u), 'measured');%SNR db scale
        
        %%%%%%%%%%%%여기까지가 기저대역 신호 생성%%%%%%%%%%%%
        cf=linspace(-fs/2, fs/2-fs/length(Rx), length(Rx));    
        
    
        for k=1:length(sample_shift)
            shifted_input=circshift(Rx, sample_shift(k)); %sample_shift는 벡터
            quadratic=Rx.*shifted_input; %x(t)x*(t-tau)
            Fm(k,:)=fftshift(fft(quadratic));
            cvm(k,:)=(1/length(quadratic))*Fm(k,:);  %covarinace vector=[real(cvm(:,i)) imag(cvm(:,i))]   
        end
        L=length(Rx)/4+1;
        w = kaiser(L,1);
        o=find(cf == 0);
        detection=zeros(1,length(cf)-L);
        %alpha=1+(L-1)/2:length(cf)-((L-1)/2) two-sided range => 그런데 양수쪽만
        %봐도 상관 없다. 
        for alpha=o+1:length(cf)-((L-1)/2)
            S = zeros(length(sample_shift), length(sample_shift));
            Scj = zeros(length(sample_shift), length(sample_shift));
            for m = 1 : length(sample_shift)
                for n = 1 : length(sample_shift)
                    fn=Fm(n,:);
                    fm=Fm(m,:);
                    fv1=flip(fn(alpha-(L-1)/2:alpha+(L-1)/2));
                    fv2=fm(alpha-(L-1)/2:alpha+(L-1)/2);
                    fv3=conj(fn(alpha-(L-1)/2:alpha+(L-1)/2));
    
                    S(m, n)=fv1.*fv2*w; % 윈도우 함수의 인덱스가 음수가 되지 않도록 조정
                    Scj(m, n)=fv2.*fv3*w;
                end
            end
            S = S / (length(Rx) * L);
            Scj = Scj / (length(Rx) * L);
            sig1=real((S+Scj)/2);
            sig2=imag((S-Scj)/2);
            sig3=imag((S+Scj)/2);
            sig4=real((Scj-S)/2);
            CM=[sig1 sig2;sig3 sig4];
            v=[real(cvm(:,alpha).') imag(cvm(:,alpha).')];
            ML=length(Rx)*v*inv(CM)*(v.');
            
            threshold=chi2inv(1-CFAR,length(v)); %자유도: Covariance vector 길이
            if ML >= threshold
                detection(alpha-(L-1)/2)=1;
            end
        end
        sensing=sensing+sign(sum(detection));
    end
    DP(u)=sensing/simulation;
end
figure(1)
plot(SNR, DP)
title("CFAR=1/(10^5), tau=0:3, Number of sample=10160")
xlabel('SNR(dB)')
ylabel('Detection Probability')