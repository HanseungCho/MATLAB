clc; clear;
M=2; N=10^4;
samplingrate_bit=100;
T_bit=1/(10^5);
Fs=samplingrate_bit/T_bit;
t=0:(1/Fs):samplingrate_bit*(N/sqrt(M))*(1/Fs)-(1/Fs);
[symbols ,Rx, Tx, noise]=FSK(M, N, T_bit, samplingrate_bit, 10, 0);
figure;
%parameter : M,N(number of bits),T_bit,samplingrate_bit,Ebn0_db,coherence
FFT_signal=abs(fft(Tx));
subplot(2,2,1);
plot(Tx(1:400)); title("Tx(1:400)");
subplot(2,2,2);
plot(Rx(1:400)); title("Rx(1:400)");
subplot(2,2,3);
F=Fs/length(FFT_signal)*(0:length(FFT_signal)-1);
plot(F,FFT_signal); title("Rx FFT");
%M&frequencies detection
threshold=max(FFT_signal)/2;
estimated_m=length(F(FFT_signal>threshold))/2;
estimated_ff=F(FFT_signal>threshold); 
estimated_f=zeros(1,estimated_m);
for p=1:estimated_m
estimated_f(p)=estimated_ff(p);
end
disp(estimated_f);
disp(append("M Value:", string(estimated_m)));

subplot(2,2,4);
figure;
BPF_Rx=zeros(estimated_m, length(Rx));
for M_num=1:estimated_m
    BPF_Rx(M_num,:)=bandpass(Rx, [estimated_f(M_num)-0.1*(estimated_f(2)-estimated_f(1)), estimated_f(M_num)+0.1*(estimated_f(2)-estimated_f(1))], Fs);
    subplot(1,estimated_m,M_num);
    plot(BPF_Rx(M_num,:));
end

%symbol estimation
[c,l] = dwt(Rx_bpf, 'haar');
Rx_dwt=idwt(c, zeros(size(c)), 'haar'); %wavelet으로 noise 완화
Rx_fft=abs(fft(Rx_dwt));
plot(Rx_fft); title("Rx DWT FFT");
%STFT demodulation
figure;
window = kaiser(100,0.5); fft_length=100;
stft(Rx_dwt(1:400), Fs, 'Window', window, 'OverlapLength', 0, 'FFTLength', fft_length);
[STFT_s, STFT_f, STFT_t] = stft(Rx_dwt, Fs, 'Window', window, 'OverlapLength', 0, 'FFTLength', fft_length);
Srow=zeros(1,length(estimated_f));
for q=1:length(estimated_f)
Srow(q)=find(STFT_f==estimated_f(q));
end
[maxS, reconsymbol] = max(STFT_s(Srow,:));
reconsymbol=reconsymbol-1; %symbol offset => 0
number_of_error = 0;
for i = 1:length(reconsymbol)
    if reconsymbol(i) ~= symbols(i)
        number_of_error = number_of_error + 1;
    end
end
disp(append("number of error: ", string(number_of_error), "   SER: ", string(number_of_error/length(symbols))));
