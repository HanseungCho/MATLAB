clc;
clear;
T=1/(10^5);
t=-10*T:T/10:10*T;%샘플링 10^6
Fs=10^6; 
b=0.5;
sym=zeros(1,1000);
data_I_QPSK=(randi([0,1],1,1000)*2-1)/sqrt(10);
data_Q_QPSK=(randi([0,1],1,1000)*2-1)/sqrt(10); % 1/4확률로 평균 
data_constant1=randi([0,1],1,1000)*2+1;
data_constant2=randi([0,1],1,1000)*2+1;
data_I=data_I_QPSK.*data_constant1;
data_Q=data_Q_QPSK.*data_constant2;

for u=1:1000
sym(u)=(data_I(u)+(1j)*data_Q(u));
end
%%%%%%%%%%%%%%%%%%oversampling
oversample_sym=zeros(1,10000);
for k=1:1000
oversample_sym(10*(k-1)+1)=sym(k);
    for l=2:10
        oversample_sym(10*(k-1)+l)=0;
    end
end

h1_t=sinc(t/T).*cos(pi*b*t/T)./(1-(2*b*t/T).^2);  
filtered=conv(oversample_sym,h1_t, 'same');
downsample_filtered=zeros(1,1000);
for m=1:1000
    downsample_filtered(m)=filtered(m*10-9);
end
%%%%%%%%%%%%%%%%%%%%%%%non offset
figure(1) %%space diagram
dI=real(downsample_filtered);
dQ=imag(downsample_filtered);
scatter(dI,dQ);
figure(2) %%eye pattern
L=length(sym);
for p=1:(L/2-1)
   R=real(filtered(p*20-19:p*20+1));
   plot(R)
   hold on
end
figure(3)
plot(filtered);
