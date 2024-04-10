clc;
clear;
N=10^6;
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
signal_d0=zeros(1,4*N);
signal_d1=zeros(1,4*N);
BER0=zeros(1,14);
SER0=zeros(1,14);
BER1=zeros(1,14);
SER1=zeros(1,14);
error_bit0=zeros(1,14);
error_sim0=zeros(1,14);
error_bit1=zeros(1,14);
error_sim1=zeros(1,14);
dI0=zeros(1,N);
dQ0=zeros(1,N);
dI1=zeros(1,N);
dQ1=zeros(1,N);
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
raych=(randn(1,N)+randn(1,N)*1j)/sqrt(2);
signal0=symbol+noise;
signal1=(symbol.*raych+noise)./raych;

for p=1:N
dI0(p)=real(signal0(p));
dQ0(p)=imag(signal0(p));
dI1(p)=real(signal1(p));
dQ1(p)=imag(signal1(p));

if (dI0(p)>(2/sqrt(10)))
    signal_d0(4*p-3)=0;
    signal_d0(4*p-1)=1;
elseif (dI0(p)<=(2/sqrt(10))) && (dI0(p)>0)
    signal_d0(4*p-3)=0;
    signal_d0(4*p-1)=0;
elseif (dI0(p)<=0) && (dI0(p)>(-2/sqrt(10)))
    signal_d0(4*p-3)=1;
    signal_d0(4*p-1)=0;
elseif (dI0(p)<=(-2/sqrt(10)))
    signal_d0(4*p-3)=1;
    signal_d0(4*p-1)=1;
end

if (dQ0(p)>(2/sqrt(10)))
    signal_d0(4*p-2)=0;
    signal_d0(4*p)=1;
elseif (dQ0(p)<=(2/sqrt(10)) && (dQ0(p)>0))
    signal_d0(4*p-2)=0;
    signal_d0(4*p)=0;
elseif (dQ0(p)<=0) && (dQ0(p)>(-2/sqrt(10)))
    signal_d0(4*p-2)=1;
    signal_d0(4*p)=0;
elseif (dQ0(p)<=(-2/sqrt(10)))
    signal_d0(4*p-2)=1;
    signal_d0(4*p)=1;
end

if (signal_d0(4*p)==symbol_bits(4*p)) && (signal_d0(4*p-1)==symbol_bits(4*p-1)) && (signal_d0(4*p-2)==symbol_bits(4*p-2)) && (signal_d0(4*p-3)==symbol_bits(4*p-3))
    error_sim0(l)=error_sim0(l)+0;
else
    error_sim0(l)=error_sim0(l)+1;
end
%%%%%%%%%%%%%%%%%%rayleigh%%%%%
if (dI1(p)>(2/sqrt(10)))
    signal_d1(4*p-3)=0;
    signal_d1(4*p-1)=1;
elseif (dI1(p)<=(2/sqrt(10))) && (dI1(p)>0)
    signal_d1(4*p-3)=0;
    signal_d1(4*p-1)=0;
elseif (dI1(p)<=0) && (dI1(p)>(-2/sqrt(10)))
    signal_d1(4*p-3)=1;
    signal_d1(4*p-1)=0;
elseif (dI1(p)<=(-2/sqrt(10)))
    signal_d1(4*p-3)=1;
    signal_d1(4*p-1)=1;
end

if (dQ1(p)>(2/sqrt(10)))
    signal_d1(4*p-2)=0;
    signal_d1(4*p)=1;
elseif (dQ1(p)<=(2/sqrt(10)) && (dQ1(p)>0))
    signal_d1(4*p-2)=0;
    signal_d1(4*p)=0;
elseif (dQ1(p)<=0) && (dQ1(p)>(-2/sqrt(10)))
    signal_d1(4*p-2)=1;
    signal_d1(4*p)=0;
elseif (dQ1(p)<=(-2/sqrt(10)))
    signal_d1(4*p-2)=1;
    signal_d1(4*p)=1;
end

if (signal_d1(4*p)==symbol_bits(4*p)) && (signal_d1(4*p-1)==symbol_bits(4*p-1)) && (signal_d1(4*p-2)==symbol_bits(4*p-2)) && (signal_d1(4*p-3)==symbol_bits(4*p-3))
    error_sim1(l)=error_sim1(l)+0;
else
    error_sim1(l)=error_sim1(l)+1;
end
end

for t=1:(4*N)
if signal_d0(t)==symbol_bits(t)
    error_bit0(l)=error_bit0(l)+0;
else 
    error_bit0(l)=error_bit0(l)+1;
end
if signal_d1(t)==symbol_bits(t)
    error_bit1(l)=error_bit1(l)+0;
else 
    error_bit1(l)=error_bit1(l)+1;
end
end

BER0(l)=error_bit0(l)/(N*4);
SER0(l)=error_sim0(l)/N;
BER1(l)=error_bit1(l)/(N*4);
SER1(l)=error_sim1(l)/N;
end
[ber0,ser0]=berawgn(ebn0_db,'qam',M,'nondiff'); 
[ber1,ser1]=berfading(ebn0_db,'qam',M,1); 
semilogy(x,BER0,'-o',x,SER0,'-s',x,BER1,'-o',x,SER1,'-s',x,ber0,x,ser0,x,ber1,x,ser1);
legend('BER_AWGN','SER_AWGN','BER_Rayleigh','SER_Rayleigh','theo_BER','theo_SER','theo_BER_Rayleigh','theo_SER_Rayleigh');
