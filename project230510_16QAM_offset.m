clc;
clear;
T=1/(10^5);
t=-10*T:T/10:10*T;%샘플링 10^6
t_offset=zeros(1,2000);
Fs=10^6; 
for e=1:1998
t_offset(e+1)=e/Fs;
end
b=0.5;
sym=zeros(1,1000);
sym_2=zeros(1,1000);
sym_3=zeros(1,1000);
sym_4=zeros(1,1000);
data_I_QPSK=(randi([0,1],1,1000)*2-1)/sqrt(10);
data_Q_QPSK=(randi([0,1],1,1000)*2-1)/sqrt(10); % 1/4확률로 평균 
data_constant1=randi([0,1],1,1000)*2+1;
data_constant2=randi([0,1],1,1000)*2+1;
data_I=data_I_QPSK.*data_constant1;
data_Q=data_Q_QPSK.*data_constant2;
for u=1:1000
sym(u)=(data_I(u)+(1j)*data_Q(u));
sym_2(u)=(data_I(u)+(1j)*data_Q(u))/sqrt(2)*exp(1j*pi/9);
sym_3(u)=(data_I(u)+(1j)*data_Q(u))/sqrt(2)*exp(1j*2*pi*1000*t_offset(u));
sym_4(u)=real(sym(u))+(1j)*imag(sym_2(u));
end
%%%%%%%%%%%%%%%%%%oversampling
oversample_sym=zeros(1,10000);
oversample_sym_2=zeros(1,10000);
oversample_sym_3=zeros(1,10000);
oversample_sym_4=zeros(1,10000);
for k=1:1000
oversample_sym(10*(k-1)+1)=sym(k);
oversample_sym_2(10*(k-1)+1)=sym_2(k);
oversample_sym_3(10*(k-1)+1)=sym_3(k);
oversample_sym_4(10*(k-1)+1)=sym_4(k);
    for l=2:10
        oversample_sym(10*(k-1)+l)=0;
        oversample_sym_2(10*(k-1)+l)=0;
        oversample_sym_3(10*(k-1)+l)=0;
        oversample_sym_4(10*(k-1)+l)=0;
    end
end

h1_t=sinc(t/T).*cos(pi*b*t/T)./(1-(2*b*t/T).^2);  
filtered=conv(oversample_sym,h1_t, 'same');
filtered_2=conv(oversample_sym_2,h1_t, 'same');
filtered_3=conv(oversample_sym_3,h1_t, 'same');
filtered_4=conv(oversample_sym_4,h1_t, 'same');
downsample_filtered=zeros(1,1000);
downsample_filtered_2=zeros(1,1000);
downsample_filtered_3=zeros(1,1000);
downsample_filtered_4=zeros(1,1000);
for m=1:1000
    downsample_filtered(m)=filtered(m*10-9);
    downsample_filtered_2(m)=filtered_2(m*10-9);
    downsample_filtered_3(m)=filtered_3(m*10-9);
    downsample_filtered_4(m)=filtered_4(m*10-9);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%phase offset = 20
figure(4) %%space diagram
dI_2=real(downsample_filtered_2);
dQ_2=imag(downsample_filtered_2);
scatter(dI_2,dQ_2);
figure(5) %%eye pattern
L=length(sym);
for p=1:(L/2-1)
   R=real(filtered_2(p*20-19:p*20+1));
   plot(R)
   hold on
end
figure(6)
plot(filtered_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%frequency offset =1000hz
figure(7) %%space diagram
dI_3=real(downsample_filtered_3);
dQ_3=imag(downsample_filtered_3);
scatter(dI_3,dQ_3);
figure(8) %%eye pattern
L=length(sym);
for p=1:(L/2-1)
   R=real(filtered_3(p*20-19:p*20+1));
   plot(R)
   hold on
end
figure(9)
plot(filtered_3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%imbalance
figure(10) %%space diagram
dI_4=real(downsample_filtered_4);
dQ_4=imag(downsample_filtered_4);
scatter(dI_4,dQ_4);
figure(11) %%eye pattern
L=length(sym);
for p=1:(L/2-1)
   R=imag(filtered_4(p*20-19:p*20+1));
   plot(R)
   hold on
end
figure(12)
plot(filtered_4);
