clc;
clear;
N=10^6;
theta=zeros(1,4);
symbol_bits=zeros(1,2*N);
ebn0_db=-3:10;
ebn0=zeros(1,14);
sig=zeros(1,14);
signal_d0=zeros(1,N);
signal_d1=zeros(1,N);
BER0=zeros(1,14);
SER0=zeros(1,14);
BER1=zeros(1,14);
SER1=zeros(1,14);
error_bit0=zeros(1,14);
error_sim0=zeros(1,14);
error_bit1=zeros(1,14);
error_sim1=zeros(1,14);
M=4;
dI0=zeros(1,N);
dQ0=zeros(1,N);
dI1=zeros(1,N);
dQ1=zeros(1,N);
for k=1:4
    theta(k)=exp(1j*(2*pi/4*(k-1)+pi/4));
end
symbol=randsrc(1,N,theta);
for u=1:N
    symbol_bits(2*u-1)=(1/2)*(1+sign(real(symbol(u))));
    symbol_bits(2*u)=(1/2)*(1+sign(imag(symbol(u))));
end
for k=1:14
ebn0(k)=10^(ebn0_db(k)/10);
sig(k)=sqrt(1/(2*2*ebn0(k)));
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

if dI0(p)>0 
    signal_d0(2*p-1)=1;
else 
    signal_d0(2*p-1)=0;
end

if dQ0(p)>0
    signal_d0(2*p)=1;
else 
    signal_d0(2*p)=0;
end

if (signal_d0(2*p-1)==symbol_bits(2*p-1)) && (signal_d0(2*p)==symbol_bits(2*p))
    error_sim0(l)=error_sim0(l)+0;
else
    error_sim0(l)=error_sim0(l)+1;
end

if dI1(p)>0 
    signal_d1(2*p-1)=1;
else 
    signal_d1(2*p-1)=0;
end

if dQ1(p)>0
    signal_d1(2*p)=1;
else 
    signal_d1(2*p)=0;
end

if (signal_d1(2*p-1)==symbol_bits(2*p-1)) && (signal_d1(2*p)==symbol_bits(2*p))
    error_sim1(l)=error_sim1(l)+0;
else
    error_sim1(l)=error_sim1(l)+1;
end

end

for t=1:(2*N)
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
%%%%%%%%%%%%%%%%%%%%%%%%%
BER0(l)=error_bit0(l)/(2*N);
SER0(l)=error_sim0(l)/N;
BER1(l)=error_bit1(l)/(2*N);
SER1(l)=error_sim1(l)/N;
end
[ber0,ser0]=berawgn(ebn0_db,'psk',M,'nondiff');
[ber1,ser1]=berfading(ebn0_db,'psk',M,1); 
semilogy(x,BER0,x,SER0,x,BER1,'-o',x,SER1,'-s',x,ber0,x,ser0,x,ber1,x,ser1);
legend('BER AWGN','SER AWGN','BER Rayleigh','SER Rayleigh','theo BER','theo SER','theo BER Rayleigh','theo SER Rayleigh');