clc;
clear;
bits=randi([0,1],1,10^6)*2-1;
ebn0_db=-3:10;
ebn0=zeros(1,14);
sig=zeros(1,14);
signal_d=zeros(1,10^6);
BER=zeros(1,14);
SER=zeros(1,14);
M=2;

for k=1:14
ebn0(k)=10^(ebn0_db(k)/10);
sig(k)=sqrt(1/(2*ebn0(k)));
end
  %sig 조정
x=1:14;
error_count=zeros(1,14);
for l=1:14
noise=randn(1,10^6)*sig(l);
signal=bits+noise;
for t=1:10^6
if signal(t)>0 
    signal_d(t)=1;
else 
    signal_d(t)=-1;
end
if signal_d(t)==bits(t)
    error_count(l)=error_count(l)+0;
else 
    error_count(l)=error_count(l)+1;
end
end
BER(l)=error_count(l)/(10^6);
SER(l)=error_count(l)/(10^6);
end
[ber,ser]=berawgn(ebn0_db,'psk',M,'nondiff'); 
semilogy(x,BER,x,SER,x,ber,x,ser);
legend('BER','SER','theo_BER','theo_SER');





