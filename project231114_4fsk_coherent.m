clc;
clear;
N=10^4; %%bits numbers
T=2/(10^5); %%symbol duration(1bit)
D=200; %한 symbol안에 sampling 개수
Fs=1/T*D;
ebn0_db=-3:10; 
ebn0=zeros(1,14); 
sig=zeros(1,14); 
BER0=zeros(1,14); 
SER0=zeros(1,14); 
error_count0=zeros(1,14);
error_count1=zeros(1,14);
M=4;
phase=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sig%%%%%%%%  
for k=1:14
ebn0(k)=10^(ebn0_db(k)/10);
sig(k)=sqrt(1/(2*ebn0(k)));
end
%%%%%%%%%%%%%%%%%%%%%TX with modulation%%%%%%%%%%%%%%%%%%%%%%
bits=randi([0,1],1,N);
t=0:(1/Fs):D*(N/sqrt(M))*(1/Fs)-(1/Fs);
t_symbol=0:(1/Fs):D*(1/Fs)-(1/Fs);
coherent_t=zeros(1,N/sqrt(M)*D);
%%%%%%%%%%%%%%%%%%%making time for coherent sinusoid%%%%%%%%%%%%%%%%%
for v=1:N/sqrt(M)
    coherent_t((v-1)* D +1 : v * D) = t(1 : D);
end
f1=Fs/50;
f2=Fs/50+1/(2*T);
f3=Fs/50+2/(2*T);
f4=Fs/50+3/(2*T);

symbol_f1=cos(2*pi*f1*t_symbol+phase);
E_f1=sum(symbol_f1.^2);
symbol_f1=symbol_f1/sqrt(E_f1/sqrt(M)); %symbol energy를 1로 맞춤.
symbol_f2=cos(2*pi*f2*t_symbol+phase);
E_f2=sum(symbol_f2.^2);
symbol_f2=symbol_f2/sqrt(E_f2/sqrt(M));
symbol_f3=cos(2*pi*f3*t_symbol+phase);
E_f3=sum(symbol_f3.^2);
symbol_f3=symbol_f3/sqrt(E_f3/sqrt(M)); 
symbol_f4=cos(2*pi*f4*t_symbol+phase);
E_f4=sum(symbol_f4.^2);
symbol_f4=symbol_f4/sqrt(E_f4/sqrt(M));

Tx = [];
for q=1:N/sqrt(M)
    if bits(2*q-1 : 2*q) == [0 0]
        Tx = [Tx symbol_f1];
    elseif bits(2*q-1 : 2*q) == [0 1]
        Tx = [Tx symbol_f2];
    elseif bits(2*q-1 : 2*q) == [1 0]
        Tx = [Tx symbol_f3];
    else 
        Tx = [Tx symbol_f4];
    end
end

SUM=zeros(14*4,N/sqrt(M));
for l=1:14
%%%%%%%%%%%%%%%%%%%%%%%%%Rx signal with demodulation%%%%%%%%%%%
noise=randn(1,length(Tx))*sig(l)+randn(1,length(Tx))*sig(l)*1j;
Rx0=Tx+noise;

Rx0_w1=real(Rx0) * sqrt(sqrt(M)/E_f1) .*cos(2*pi*f1*coherent_t+phase);
Rx0_w2=real(Rx0) * sqrt(sqrt(M)/E_f2) .*cos(2*pi*f2*coherent_t+phase);
Rx0_w3=real(Rx0) * sqrt(sqrt(M)/E_f3) .*cos(2*pi*f3*coherent_t+phase);
Rx0_w4=real(Rx0) * sqrt(sqrt(M)/E_f4) .*cos(2*pi*f4*coherent_t+phase);
Rx0_bits=zeros(1,N);
for k=1:N/sqrt(M)
    Rx0_w1_symbol=Rx0_w1((k-1)*D+1 : k*D);
    Rx0_w2_symbol=Rx0_w2((k-1)*D+1 : k*D);
    Rx0_w3_symbol=Rx0_w3((k-1)*D+1 : k*D);
    Rx0_w4_symbol=Rx0_w4((k-1)*D+1 : k*D);
    SUM(4*(l-1)+1,k)=sum(Rx0_w1_symbol);
    SUM(4*(l-1)+2,k)=sum(Rx0_w2_symbol);
    SUM(4*(l-1)+3,k)=sum(Rx0_w3_symbol);
    SUM(4*(l-1)+4,k)=sum(Rx0_w4_symbol);
    E1=sum(Rx0_w1_symbol); E3=sum(Rx0_w3_symbol);
    E2=sum(Rx0_w2_symbol); E4=sum(Rx0_w4_symbol);
    %%%%%%%%%%%%AWGN decision%%%%%%%%%%%%%%%%%%
    if max([E1, E2, E3, E4]) == E1
        Rx0_bits(2*k-1 : 2*k) = [0 0];
    elseif max([E1, E2, E3, E4]) == E2
        Rx0_bits(2*k-1 : 2*k) = [0 1];
    elseif max([E1, E2, E3, E4]) == E3
        Rx0_bits(2*k-1 : 2*k) = [1 0];
    else
        Rx0_bits(2*k-1 : 2*k) = [1 1];
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
[ber0,ser0]=berawgn(ebn0_db,'fsk',M,'coherent'); 
semilogy(x,BER0,x,SER0,x,ber0,'-o',x,ser0,'-s');
legend('BER AWGN','SER AWGN','theo AWGN BER','theo AWGN SER');



