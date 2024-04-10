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
error_count0=zeros(1,14);
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
f1=Fs/100;
f2=Fs/100+1/(T);
phase=pi/3;
freq=zeros(1,N*D);
for q=1:N
    if bits(q) == 0
        freq((q-1)*D +1 : q*D) = f1;
    else
        freq((q-1)*D +1 : q*D) = f2;
    end
end
Tx_non_normalized=cos(2*pi*freq.*t+phase);
Inphase_f1       =cos(2*pi*f1*t);
Quadrature_f1    =sin(2*pi*f1*t);
Inphase_f2       =cos(2*pi*f2*t);
Quadrature_f2    =sin(2*pi*f2*t);
for d=1:N
    symbol_E              =sum(Tx_non_normalized((d-1)*D+1 : d*D).^2);
    symbol_E_Inphase_f1   =sum(Inphase_f1((d-1)*D+1 : d*D).^2);
    symbol_E_Quadrature_f1=sum(Quadrature_f1((d-1)*D+1 : d*D).^2);
    symbol_E_Inphase_f2   =sum(Inphase_f2((d-1)*D+1 : d*D).^2);
    symbol_E_Quadrature_f2=sum(Quadrature_f2((d-1)*D+1 : d*D).^2);
    Tx((d-1)*D+1 : d*D) = Tx_non_normalized((d-1)*D+1 : d*D)/sqrt(symbol_E);
    Inphase_f1((d-1)*D+1 : d*D) = Inphase_f1((d-1)*D+1 : d*D)/sqrt(symbol_E_Inphase_f1);
    Quadrature_f1((d-1)*D+1 : d*D) = Quadrature_f1((d-1)*D+1 : d*D)/sqrt(symbol_E_Quadrature_f1);
    Inphase_f2((d-1)*D+1 : d*D) = Inphase_f2((d-1)*D+1 : d*D)/sqrt(symbol_E_Inphase_f2);
    Quadrature_f2((d-1)*D+1 : d*D) = Quadrature_f2((d-1)*D+1 : d*D)/sqrt(symbol_E_Quadrature_f2);
end
for l=1:14
%%%%%%%%%%%%%%%%%%%%%%%%%Rx signal with demodulation%%%%%%%%%%%
noise=randn(1,length(Tx))*sig(l)+1j*randn(1,length(Tx))*sig(l);
Rx0=Tx+noise;

Rx0_w1I=real(Rx0) .*Inphase_f1;  %이 부분은 에너지 노멀라이즈 안해도 결과는 동일함.
Rx0_w1Q=real(Rx0) .*Quadrature_f1;
Rx0_w2I=real(Rx0) .*Inphase_f2;
Rx0_w2Q=real(Rx0) .*Quadrature_f2;
Rx0_bits=zeros(1,N);

for k=1:N
    Rx0_w1I_symbol=Rx0_w1I((k-1)*D+1 : k*D);
    Rx0_w1Q_symbol=Rx0_w1Q((k-1)*D+1 : k*D);
    Rx0_w2I_symbol=Rx0_w2I((k-1)*D+1 : k*D);
    Rx0_w2Q_symbol=Rx0_w2Q((k-1)*D+1 : k*D);
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
end
BER0(l)=error_count0(l)/N;
end

x=-3:10;
[ber0,ser0]=berawgn(ebn0_db,'fsk',M,'noncoherent'); 
semilogy(x,BER0,x,ber0,'-o');
legend('BER AWGN','theo AWGN BER');



