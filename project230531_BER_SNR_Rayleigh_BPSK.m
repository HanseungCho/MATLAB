clc;
clear;
N=10^6;
bits=randi([0,1],1,N)*2-1;
ebn0_db=-3:10;
ebn0=zeros(1,14);
sig=zeros(1,14);
signal_d0=zeros(1,N);
signal_d1=zeros(1,N);
BER0=zeros(1,14);
SER0=zeros(1,14);
BER1=zeros(1,14);
SER1=zeros(1,14);
M=2;

for k=1:14
ebn0(k)=10^(ebn0_db(k)/10);
sig(k)=sqrt(1/(2*ebn0(k)));
end

%sig 조정
x=1:14;
error_count0=zeros(1,14);
error_count1=zeros(1,14);
for l=1:14
noise=randn(1,N)*sig(l)+1j*randn(1,N)*sig(l);
raych=(randn(1,N)+randn(1,N)*1j)/sqrt(2);
signal0=bits+noise;
signal1=(bits.*raych+noise)./raych;
for t=1:N
if signal0(t)>0 
    signal_d0(t)=1;
else 
    signal_d0(t)=-1;
end
if signal_d0(t)==bits(t)
    error_count0(l)=error_count0(l)+0;
else 
    error_count0(l)=error_count0(l)+1;
end
end
for t=1:N
if signal1(t)>0 
    signal_d1(t)=1;
else 
    signal_d1(t)=-1;
end
if signal_d1(t)==bits(t)
    error_count1(l)=error_count1(l)+0;
else 
    error_count1(l)=error_count1(l)+1;
end
end
BER0(l)=error_count0(l)/N;
SER0(l)=error_count0(l)/N;
BER1(l)=error_count1(l)/N;
SER1(l)=error_count1(l)/N;
end
[ber0,ser0]=berawgn(ebn0_db,'psk',M,'nondiff'); 
[ber1,ser1]=berfading(ebn0_db,'psk',M,1); 
semilogy(x,BER0,'-o',x,BER1,'-s',x,ber0,x,ber1);
legend('BER_AWGN','BER_Rayleigh','theo_BER','theo_SER');





