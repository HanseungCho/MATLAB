clc;
clear;
theta=zeros(1,4);
symbol_bits=zeros(1,10^6);
ebn0_db=-3:10;
ebn0=zeros(1,14);
sig=zeros(1,14);
signal_d=zeros(1,10^6);
BER=zeros(1,14);
SER=zeros(1,14);
error_bit=zeros(1,14);
error_sim=zeros(1,14);
M=4;
dI=zeros(1,(10^6)/2);
dQ=zeros(1,(10^6)/2);
for k=1:4
    theta(k)=exp(1j*(2*pi/4*(k-1)+pi/4));
end
symbol=randsrc(1,(10^6)/2,theta);
for u=1:((10^6)/2)
    symbol_bits(2*u-1)=(1/2)*(1+sign(real(symbol(u))));
    symbol_bits(2*u)=(1/2)*(1+sign(imag(symbol(u))));
end
for k=1:14
ebn0(k)=10^(ebn0_db(k)/10);
sig(k)=sqrt(1/(2*2*ebn0(k)));
end
x=1:14;

for l=1:14
noise=randn(1,(10^6)/2)*sig(l)+randn(1,(10^6)/2)*sig(l)*1j;
signal=symbol+noise;
for p=1:((10^6)/2)
dI(p)=real(signal(p));
dQ(p)=imag(signal(p));

if dI(p)>0 
    signal_d(2*p-1)=1;
else 
    signal_d(2*p-1)=0;
end

if dQ(p)>0
    signal_d(2*p)=1;
else 
    signal_d(2*p)=0;
end

if (signal_d(2*p-1)==symbol_bits(2*p-1)) && (signal_d(2*p)==symbol_bits(2*p))
    error_sim(l)=error_sim(l)+0;
else
    error_sim(l)=error_sim(l)+1;
end

end

for t=1:(10^6)
if signal_d(t)==symbol_bits(t)
    error_bit(l)=error_bit(l)+0;
else 
    error_bit(l)=error_bit(l)+1;
end
end

BER(l)=error_bit(l)/(10^6);
SER(l)=error_sim(l)/((10^6)/2);
end
[ber,ser]=berawgn(ebn0_db,'psk',M,'nondiff'); 
semilogy(x,BER,x,SER,x,ber,x,ser);
legend('BER','SER','theo_BER','theo_SER');
