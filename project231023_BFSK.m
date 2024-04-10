clc;
clear;
N=10^5; %%bits numbers
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
t=0:(1/Fs):1-(1/Fs);
coherent_t=zeros(1,N*D);
%%%%%%%%%%%%%%%%%%%making time for coherent sinusoid%%%%%%%%%%%%%%%%%
for v=1:N
    coherent_t((v-1)* D +1 : v * D) = t(1 : D);
end
f1=1*10^5;
f2=3*10^5;
freq=zeros(1,N*10);
for q=1:N
    if bits(q) == 0
        freq((q-1)*D +1 : q*D) = f1;
    else
        freq((q-1)*D +1 : q*D) = f2;
    end
end
Tx=sqrt(2/D)*cos(2*pi*freq.*coherent_t);
SUM=zeros(14*4,N);
for l=1:14
%%%%%%%%%%%%%%%%%%%%%%%%%Rx signal with demodulation%%%%%%%%%%%
noise=sqrt(1/2)*(randn(1,length(Tx))*sig(l)+randn(1,length(Tx))*sig(l));
raych=sqrt(1/2)*(randn(1,length(Tx))+randn(1,length(Tx))*1j);
Rx0=Tx+noise;
Rx1=(Tx.*raych+noise)./raych;

Rx0_w1=real(Rx0) * sqrt(2/D) .*cos(2*pi*f1*coherent_t);
Rx0_w2=real(Rx0) * sqrt(2/D) .*cos(2*pi*f2*coherent_t);
Rx1_w1=real(Rx1) * sqrt(2/D) .*cos(2*pi*f1*coherent_t);
Rx1_w2=real(Rx1) * sqrt(2/D) .*cos(2*pi*f2*coherent_t);
Rx0_bits=zeros(1,N);
Rx1_bits=zeros(1,N);
for k=1:N
    Rx0_w1_symbol=Rx0_w1((k-1)*D+1 : k*D);
    Rx0_w2_symbol=Rx0_w2((k-1)*D+1 : k*D);
    Rx1_w1_symbol=Rx1_w1((k-1)*D+1 : k*D);
    Rx1_w2_symbol=Rx1_w2((k-1)*D+1 : k*D);
    SUM(4*(l-1)+1,k)=sum(Rx0_w1_symbol);
    SUM(4*(l-1)+2,k)=sum(Rx0_w2_symbol);
    SUM(4*(l-1)+3,k)=sum(Rx1_w1_symbol);
    SUM(4*(l-1)+4,k)=sum(Rx1_w2_symbol);
    %%%%%%%%%%%%AWGN decision%%%%%%%%%%%%%%%%%%
    if sum(Rx0_w1_symbol) > sum(Rx0_w2_symbol)
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
    if sum(Rx1_w1_symbol) > sum(Rx1_w2_symbol)
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
[ber0,ser0]=berawgn(ebn0_db,'fsk',M,'coherent'); 
[ber1,ser1]=berfading(ebn0_db,'fsk',M,1,'coherent'); 
semilogy(x,BER0,'-o',x,BER1,x,ber0,'-o',x,ber1,'-s');
legend('BER AWGN','BER Rayleigh','theo AWGN BER','theo Rayleigh BER');



