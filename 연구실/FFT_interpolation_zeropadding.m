clc; clear;
N=10^3;
Tb=1/(10^3);
fs=100/(Tb); %oversampling 100배
bits=randi([0,1], 1, N);
for k=1:length(bits)
    oversampled_bits((k-1)*100+1:k*100)=bits(k);
end
NRZ=oversampled_bits*2-1;
Bipolar_RZ=NRZ;
for k=1:length(bits)
    Bipolar_RZ((k-1)*100+51:k*100)=0;
end
Bi_Phase_L_unit=[ones(1,50), -1*ones(1,50)];
for k=1:1000
    if bits(k)==1
        Bi_Phase_L((k-1)*100+1:k*100)=Bi_Phase_L_unit;
    else 
        Bi_Phase_L((k-1)*100+1:k*100)=-1*Bi_Phase_L_unit;
    end
end
t=0:Tb/100:Tb*10-(Tb/100);
figure(1);
subplot(4,1,1);
stem(bits(1:11));
title("message bits")
subplot(4,1,2);
plot(t,NRZ(1:10*100));
title("NRZ")
subplot(4,1,3)
plot(t,Bipolar_RZ(1:10*100));
title("Bipolar RZ")
subplot(4,1,4)
plot(t,Bi_Phase_L(1:10*100));
title("Bi phase L")
figure(2);
NRZ_fft=abs(fftshift(fft(NRZ)));
Bipolar_RZ_fft=abs(fftshift(fft(Bipolar_RZ)));
Bi_Phase_L_fft=abs(fftshift(fft(Bi_Phase_L)));
hold on
f=linspace(-fs/2,fs/2,N*100);
plot(f, NRZ_fft);
plot(f, Bipolar_RZ_fft);
plot(f, Bi_Phase_L_fft);


