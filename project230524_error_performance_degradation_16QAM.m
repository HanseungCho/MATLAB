clc;
clear;
N=10^5;
M=16;
w=4; %시그마에 k값 
sym=zeros(1,N);
data_I_QPSK=(randi([0,1],1,N)*2-1)/sqrt(10);
data_Q_QPSK=(randi([0,1],1,N)*2-1)/sqrt(10); % 1/4확률로 평균 
data_constant1=randi([0,1],1,N)*2+1;
data_constant2=randi([0,1],1,N)*2+1;
data_I=data_I_QPSK.*data_constant1;
data_Q=data_Q_QPSK.*data_constant2;
for u=1:N
sym(u)=(data_I(u)+(1j)*data_Q(u));
end

theta=zeros(1,4);
symbol_bits=zeros(1,4*N);
ebn0_db=-3:10;
ebn0=zeros(1,14);
sig=zeros(1,14);
signal_d=zeros(1,4*N);
BER=zeros(1,14);
SER=zeros(1,14);
error_bit=zeros(1,14);
error_sim=zeros(1,14);
dI=zeros(1,N);
dQ=zeros(1,N);
symbol=sym;

for u=1:(N)
    if (real(symbol(u)) == (3/sqrt(10)))
    symbol_bits(4*u-3)=0;
    symbol_bits(4*u-1)=1;
    elseif (real(symbol(u)) == (1/sqrt(10)))
    symbol_bits(4*u-3)=0;
    symbol_bits(4*u-1)=0;
    elseif (real(symbol(u)) == (-1/sqrt(10)))
    symbol_bits(4*u-3)=1;
    symbol_bits(4*u-1)=0;
    elseif (real(symbol(u)) == (-3/sqrt(10)))
    symbol_bits(4*u-3)=1;
    symbol_bits(4*u-1)=1;
    end

    if (imag(symbol(u)) == (3/sqrt(10)))
    symbol_bits(4*u-2)=0;
    symbol_bits(4*u)=1;
    elseif (imag(symbol(u)) == (1/sqrt(10)))
    symbol_bits(4*u-2)=0;
    symbol_bits(4*u)=0;
    elseif (imag(symbol(u)) == (-1/sqrt(10)))
    symbol_bits(4*u-2)=1;
    symbol_bits(4*u)=0;
    elseif (imag(symbol(u)) == (-3/sqrt(10)))
    symbol_bits(4*u-2)=1;
    symbol_bits(4*u)=1;
    end
end
for k=1:14
ebn0(k)=10^(ebn0_db(k)/10);
sig(k)=sqrt(1/(2*w*ebn0(k)));
end
x=1:14;

for l=1:14
noise=randn(1,N)*sig(l)+randn(1,N)*sig(l)*1j;
signal=symbol+noise;
for p=1:N
dI(p)=real(signal(p));
dQ(p)=imag(signal(p));

if (dI(p)>(2/sqrt(10)))
    signal_d(4*p-3)=0;
    signal_d(4*p-1)=1;
elseif (dI(p)<=(2/sqrt(10))) && (dI(p)>0)
    signal_d(4*p-3)=0;
    signal_d(4*p-1)=0;
elseif (dI(p)<=0) && (dI(p)>(-2/sqrt(10)))
    signal_d(4*p-3)=1;
    signal_d(4*p-1)=0;
elseif (dI(p)<=(-2/sqrt(10)))
    signal_d(4*p-3)=1;
    signal_d(4*p-1)=1;
end

if (dQ(p)>(2/sqrt(10)))
    signal_d(4*p-2)=0;
    signal_d(4*p)=1;
elseif (dQ(p)<=(2/sqrt(10)) && (dQ(p)>0))
    signal_d(4*p-2)=0;
    signal_d(4*p)=0;
elseif (dQ(p)<=0) && (dQ(p)>(-2/sqrt(10)))
    signal_d(4*p-2)=1;
    signal_d(4*p)=0;
elseif (dQ(p)<=(-2/sqrt(10)))
    signal_d(4*p-2)=1;
    signal_d(4*p)=1;
end

if (signal_d(4*p)==symbol_bits(4*p)) && (signal_d(4*p-1)==symbol_bits(4*p-1)) && (signal_d(4*p-2)==symbol_bits(4*p-2)) && (signal_d(4*p-3)==symbol_bits(4*p-3))
    error_sim(l)=error_sim(l)+0;
else
    error_sim(l)=error_sim(l)+1;
end

end

for t=1:(4*N)
if signal_d(t)==symbol_bits(t)
    error_bit(l)=error_bit(l)+0;
else 
    error_bit(l)=error_bit(l)+1;
end
end

BER(l)=error_bit(l)/(N*4);
SER(l)=error_sim(l)/N;
end
[ber,ser]=berawgn(ebn0_db,'qam',M,'nondiff'); 
semilogy(x,BER,x,SER,x,ber,x,ser);
legend('BER','SER','theo_BER','theo_SER');
