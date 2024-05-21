clc;
clear;
N=10^4; %%bits numbers
T=1/(10^5); %%symbol duration(1bit)
D=100; %한 symbol안에 sampling 개수
Fs=1/T*D;
ebn0_db=-3:10; 
ebn0=zeros(1,14);
sig=zeros(1,14);
BER0=zeros(1,14);
SER0=zeros(1,14);
BER1=zeros(1,14);
SER1=zeros(1,14);
error_count0=zeros(1,14);
error_count1=zeros(1,14);
M=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sig%%%%%%%%
for k=1:14
ebn0(k)=10^(ebn0_db(k)/10);
sig(k)=sqrt(1/(2*ebn0(k)));
end
%%%%%%%%%%%%%%%%%%%%%TX with modulation%%%%%%%%%%%%%%%%%%%%%%
bits=randi([0,1],1,N);
t=0:(1/Fs):D*N*(1/Fs)-(1/Fs);
%%%%%%%%%%%%%%%%%%%making time for coherent sinusoid%%%%%%%%%%%%%%%%%
f1=Fs/50;
f2=Fs/50+10/(T);
phase= pi/2;
freq=zeros(1,N*D);
for q=1:N
    if bits(q) == 0
        freq((q-1)*D +1 : q*D) = f1;
    else
        freq((q-1)*D +1 : q*D) = f2;
    end
end
Tx=sqrt(2/D)*cos(2*pi*freq.*t+phase);
for l=1:14
%%%%%%%%%%%%%%%%%%%%%%%%%Rx signal with demodulation%%%%%%%%%%%
noise=randn(1,length(Tx))*sig(l);
raych=sqrt(1/2)*(randn(1,length(Tx))+randn(1,length(Tx))*1j);
Rx0=Tx+noise;
Rx1=(Tx.*raych+noise)./raych;

Rx0_w1I=real(Rx0) * sqrt(2/D) .*cos(2*pi*f1*t);
Rx0_w1Q=real(Rx0) * sqrt(2/D) .*sin(2*pi*f1*t);
Rx0_w2I=real(Rx0) * sqrt(2/D) .*cos(2*pi*f2*t);
Rx0_w2Q=real(Rx0) * sqrt(2/D) .*sin(2*pi*f2*t);
Rx1_w1I=real(Rx1) * sqrt(2/D) .*cos(2*pi*f1*t);
Rx1_w1Q=real(Rx1) * sqrt(2/D) .*sin(2*pi*f1*t);
Rx1_w2I=real(Rx1) * sqrt(2/D) .*cos(2*pi*f2*t);
Rx1_w2Q=real(Rx1) * sqrt(2/D) .*sin(2*pi*f2*t);

Rx0_bits=zeros(1,N);
Rx1_bits=zeros(1,N);
for k=1:N
    Rx0_w1I_symbol=Rx0_w1I((k-1)*D+1 : k*D);
    Rx0_w1Q_symbol=Rx0_w1Q((k-1)*D+1 : k*D);
    Rx0_w2I_symbol=Rx0_w2I((k-1)*D+1 : k*D);
    Rx0_w2Q_symbol=Rx0_w2Q((k-1)*D+1 : k*D);
    Rx1_w1I_symbol=Rx1_w1I((k-1)*D+1 : k*D);
    Rx1_w1Q_symbol=Rx1_w1Q((k-1)*D+1 : k*D);
    Rx1_w2I_symbol=Rx1_w2I((k-1)*D+1 : k*D);
    Rx1_w2Q_symbol=Rx1_w2Q((k-1)*D+1 : k*D);
    %%%%%%%%%%%%AWGN decision%%%%%%%%%%%%%%%%%%
    if sum(Rx0_w1I_symbol)^2+sum(Rx0_w1Q_symbol)^2 > sum(Rx0_w2I_symbol)^2+sum(Rx0_w2Q_symbol)^2
        Rx0_bits(k)=0;
    else 
        Rx0_bits(k)=1;
    end
    if Rx0_bits(k)==bits(k)
        error_count0(l)=error_count0(l)+0;
    else 
        error_count0(l)=error_count0(l)+1;
    end
    %%%%%%%%%%%%Rayleigh decision%%%%%%%%%%%%%%%%%%
    if sum(Rx1_w1I_symbol)^2+sum(Rx1_w1Q_symbol)^2 > sum(Rx1_w2I_symbol)^2+sum(Rx1_w2Q_symbol)^2
        Rx1_bits(k)=0;
    else 
        Rx1_bits(k)=1;
    end
    if Rx1_bits(k)==bits(k)
        error_count1(l)=error_count1(l)+0;
    else 
        error_count1(l)=error_count1(l)+1;
    end
end
BER0(l)=error_count0(l)/N;
%SER0(l)=error_count0(l)/N;
BER1(l)=error_count1(l)/N;
%SER1(l)=error_count1(l)/N;
end

x=-3:10;
[ber0,ser0]=berawgn(ebn0_db,'fsk',M,'noncoherent'); 
[ber1,ser1]=berfading(ebn0_db,'fsk',M,1,'noncoherent'); 
semilogy(x,BER0,'-o',x,BER1,x,ber0,'-o',x,ber1,'-s');
legend('BER AWGN','BER Rayleigh','theo AWGN BER','theo Rayleigh BER');



