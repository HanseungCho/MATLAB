clc;
%% 
clear;
T=1/(10^5);
t=-10*T:T/10:10*T; %샘플링 10^6
Fs=10^6; 
b=0.5;
M=4;
%%M개의 경우(phase) 중 random한 복소평면dI,dQ조합 뽑기
theta=zeros(1,M);
for k=1:M
    theta(k)=exp(1j*(2*pi/M*(k-1)+pi/M));
end
data=randsrc(1,1000,theta);

%%oversampling
oversample_data=zeros(1,10000);
for k=1:1000
oversample_data(10*(k-1)+1)=data(k);
    for l=2:10
        oversample_data(10*(k-1)+l)=0;
    end
end

h1_t=sinc(t/T).*cos(pi*b*t/T)./(1-(2*b*t/T).^2);  
filtered=conv(oversample_data,h1_t, 'same');
downsample_filtered=zeros(1,1000);
for m=1:1000
    downsample_filtered(m)=filtered(m*10-9);
end

figure(1) %%eye pattern
dI=real(downsample_filtered);
dQ=imag(downsample_filtered);
scatter(dI,dQ);

figure(2) %%space diagram
L=length(data);
for p=1:(L/2-1)
   R=real(filtered(p*20-19:p*20+1));
   plot(R)
   hold on
end

figure(3)
plot(filtered); %%4번





