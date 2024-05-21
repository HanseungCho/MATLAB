clc;
clear;
N=10^5; %%bits numbers
T=2/(10^5); %%symbol duration(2bit)
D=200; %한 symbol안에 sampling 개수
Fs=1/T*D;
ebn0_db=-3:10; 
ebn0=zeros(1,14);
sig=zeros(1,14);
BER0=zeros(1,14);
SER0=zeros(1,14);
SER1=zeros(1,14);
error_count0=zeros(1,14);
error_count1=zeros(1,14);
error_count2=zeros(1,14);
M=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sig%%%%%%%%
for k=1:14
ebn0(k)=10^(ebn0_db(k)/10); 
sig(k)=sqrt(1/(2*ebn0(k)));
end
%%%%%%%%%%%%%%%%%%%%%TX with modulation%%%%%%%%%%%%%%%%%%%%%%
bits=randi([0,1],1,N);
t=0:(1/Fs):D*(N/sqrt(M))*(1/Fs)-(1/Fs);
%%%%%%%%%%%%%%%%%%%making time for coherent sinusoid%%%%%%%%%%%%%%%%%
f1=Fs/100;
f2=Fs/100+1/(T);
f3=Fs/100+2/(T);
f4=Fs/100+3/(T);
phase=pi/3;
freq=zeros(1,N/sqrt(M)*D);
for q=1:N/sqrt(M)
    if bits(2*q-1 : 2*q) == [0 0]
        freq((q-1)*D +1 : q*D) = f1;
    elseif bits(2*q-1 : 2*q) == [0 1]
        freq((q-1)*D +1 : q*D) = f2;
    elseif bits(2*q-1 : 2*q) == [1 0]
        freq((q-1)*D +1 : q*D) = f3;
    else
        freq((q-1)*D +1 : q*D) = f4;
    end
end

Tx_non_normalized=cos(2*pi*freq.*t+phase);
Inphase_f1       =cos(2*pi*f1*t);  Inphase_f3       =cos(2*pi*f3*t);
Quadrature_f1    =sin(2*pi*f1*t);  Quadrature_f3    =sin(2*pi*f3*t);
Inphase_f2       =cos(2*pi*f2*t);  Inphase_f4       =cos(2*pi*f4*t);
Quadrature_f2    =sin(2*pi*f2*t);  Quadrature_f4    =sin(2*pi*f4*t);
for d=1:N/sqrt(M)
    symbol_E              =sum(Tx_non_normalized((d-1)*D+1 : d*D).^2); %2*Eb
    symbol_E_Inphase_f1   =sum(Inphase_f1((d-1)*D+1 : d*D).^2);
    symbol_E_Quadrature_f1=sum(Quadrature_f1((d-1)*D+1 : d*D).^2);
    symbol_E_Inphase_f2   =sum(Inphase_f2((d-1)*D+1 : d*D).^2);
    symbol_E_Quadrature_f2=sum(Quadrature_f2((d-1)*D+1 : d*D).^2);
    symbol_E_Inphase_f3   =sum(Inphase_f3((d-1)*D+1 : d*D).^2);
    symbol_E_Quadrature_f3=sum(Quadrature_f3((d-1)*D+1 : d*D).^2);
    symbol_E_Inphase_f4   =sum(Inphase_f4((d-1)*D+1 : d*D).^2);
    symbol_E_Quadrature_f4=sum(Quadrature_f4((d-1)*D+1 : d*D).^2);

    Tx((d-1)*D+1 : d*D) = Tx_non_normalized((d-1)*D+1 : d*D)/sqrt(symbol_E/sqrt(M));
    Inphase_f1((d-1)*D+1 : d*D) = Inphase_f1((d-1)*D+1 : d*D)/sqrt(symbol_E_Inphase_f1/sqrt(M));
    Quadrature_f1((d-1)*D+1 : d*D) = Quadrature_f1((d-1)*D+1 : d*D)/sqrt(symbol_E_Quadrature_f1/sqrt(M));
    Inphase_f2((d-1)*D+1 : d*D) = Inphase_f2((d-1)*D+1 : d*D)/sqrt(symbol_E_Inphase_f2/sqrt(M));
    Quadrature_f2((d-1)*D+1 : d*D) = Quadrature_f2((d-1)*D+1 : d*D)/sqrt(symbol_E_Quadrature_f2/sqrt(M));
    Inphase_f3((d-1)*D+1 : d*D) = Inphase_f3((d-1)*D+1 : d*D)/sqrt(symbol_E_Inphase_f3/sqrt(M));
    Quadrature_f3((d-1)*D+1 : d*D) = Quadrature_f3((d-1)*D+1 : d*D)/sqrt(symbol_E_Quadrature_f3/sqrt(M));
    Inphase_f4((d-1)*D+1 : d*D) = Inphase_f4((d-1)*D+1 : d*D)/sqrt(symbol_E_Inphase_f4/sqrt(M));
    Quadrature_f4((d-1)*D+1 : d*D) = Quadrature_f4((d-1)*D+1 : d*D)/sqrt(symbol_E_Quadrature_f4/sqrt(M));
end
for l=1:14
%%%%%%%%%%%%%%%%%%%%%%%%%Rx signal with demodulation%%%%%%%%%%%
noise=randn(1,length(Tx))*sig(l)+1j*randn(1,length(Tx))*sig(l);
Rx0=Tx+noise;

Rx0_w1I=real(Rx0) .*Inphase_f1;  %이 부분은 에너지 노멀라이즈 안해도 결과는 동일함.
Rx0_w1Q=real(Rx0) .*Quadrature_f1;
Rx0_w2I=real(Rx0) .*Inphase_f2;
Rx0_w2Q=real(Rx0) .*Quadrature_f2;
Rx0_w3I=real(Rx0) .*Inphase_f3;  
Rx0_w3Q=real(Rx0) .*Quadrature_f3;
Rx0_w4I=real(Rx0) .*Inphase_f4;
Rx0_w4Q=real(Rx0) .*Quadrature_f4;
Rx0_bits=zeros(1,N);

for k=1:N/sqrt(M)
    Rx0_w1I_symbol=Rx0_w1I((k-1)*D+1 : k*D);
    Rx0_w1Q_symbol=Rx0_w1Q((k-1)*D+1 : k*D);
    Rx0_w2I_symbol=Rx0_w2I((k-1)*D+1 : k*D);
    Rx0_w2Q_symbol=Rx0_w2Q((k-1)*D+1 : k*D);
    Rx0_w3I_symbol=Rx0_w3I((k-1)*D+1 : k*D);
    Rx0_w3Q_symbol=Rx0_w3Q((k-1)*D+1 : k*D);
    Rx0_w4I_symbol=Rx0_w4I((k-1)*D+1 : k*D);
    Rx0_w4Q_symbol=Rx0_w4Q((k-1)*D+1 : k*D);
    %%%%%%%%%%%%AWGN decision%%%%%%%%%%%%%%%%%%
    E1=sum(Rx0_w1I_symbol)^2+sum(Rx0_w1Q_symbol)^2;
    E2=sum(Rx0_w2I_symbol)^2+sum(Rx0_w2Q_symbol)^2;
    E3=sum(Rx0_w3I_symbol)^2+sum(Rx0_w3Q_symbol)^2;
    E4=sum(Rx0_w4I_symbol)^2+sum(Rx0_w4Q_symbol)^2;
    if max([E1, E2, E3, E4]) == E1
        Rx0_bits(2*k-1 : 2*k)= [0 0];
    elseif max([E1, E2, E3, E4]) == E2
        Rx0_bits(2*k-1 : 2*k)= [0 1];
    elseif max([E1, E2, E3, E4]) == E3
        Rx0_bits(2*k-1 : 2*k)= [1 0];
    else
        Rx0_bits(2*k-1 : 2*k)= [1 1];
    end

    if Rx0_bits(2*k-1 : 2*k) == bits(2*k-1 : 2*k)
        error_count1(l)=error_count1(l)+0;
    else
        error_count1(l)=error_count1(l)+1;
    end
end
error_count0(l)=sum((Rx0_bits~=bits));
BER0(l)=error_count0(l)/N;
SER0(l)=error_count1(l)/(N/sqrt(M));
end

x=-3:10;
[ber0,ser0]=berawgn(ebn0_db,'fsk',M,'noncoherent'); 
semilogy(x,BER0,x,SER0,x,ber0,'-o',x,ser0,'-s');
legend('BER AWGN','SER AWGN','theo AWGN BER','theo AWGN SER');



