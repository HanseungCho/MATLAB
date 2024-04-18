clc;
clear;
M=4;
N=1*10^4; %%bits numbers
bits_for_symbol=log2(M);
T=bits_for_symbol/(10^5); %%symbol duration(1bit)
D=bits_for_symbol*100; %한 symbol안에 sampling 개수
Fs=1/T*D;
%f=[Fs/100+2/(T) Fs/100+3/(T)];%2FSK의 주파수 결정
f=[Fs/100 Fs/100+1/(T) Fs/100+2/(T) Fs/100+3/(T)];%4FSK의 주파수 결정
%f=[Fs/100 Fs/100+1/(T) Fs/100+2/(T) Fs/100+3/(T) Fs/100+4/(T) Fs/100+5/(T) Fs/100+6/(T) Fs/100+7/(T)];%8FSK의 주파수 결정
ebn0_db=-3:10; 
esn0=zeros(1,14);
sig=zeros(1,14);
BER0=zeros(1,14);
SER0=zeros(1,14);
bit_error_count=zeros(1,14);
sym_error_count=zeros(1,14);

%심볼 단위로 에너지를 fix해주고 노이즈를 변화시키므로 Ebn0를 Esn0로 변환시켜주어야 함.
for k=1:14
    esn0_db=ebn0_db+10*log10(bits_for_symbol); 
    esn0(k)=10^(esn0_db(k)/10); 
    sig(k)=sqrt(1/(2*esn0(k)));
end
%%%%%%%%%%%%%%%%%%%%%TX with modulation%%%%%%%%%%%%%%%%%%%%%%
bits=randi([0,1],1,N);
symbols_bits = reshape(bits, bits_for_symbol, []).';
symbols_decimal=bi2de(symbols_bits, 'left-msb');
t=0:(1/Fs):D*(N/bits_for_symbol)*(1/Fs)-(1/Fs);
%%%%%%%%%%%%%%%%%%%making time for coherent sinusoid%%%%%%%%%%%%%%%%%
phase=pi/3;
freq=zeros(1,N/bits_for_symbol*D);
for q=1:length(symbols_decimal)
        freq((q-1)*D +1 : q*D) = f(symbols_decimal(q)+1);
end
Tx_non_normalized=cos(2*pi*freq.*t+phase);
Inphase_f=zeros(M,length(Tx_non_normalized));

for d=1:N/bits_for_symbol
    symbol_E=sum(Tx_non_normalized((d-1)*D+1 : d*D).^2);
    Tx((d-1)*D+1 : d*D) = Tx_non_normalized((d-1)*D+1 : d*D)/sqrt(symbol_E);
end

for g=1:length(f)
    Inphase_f_mem=cos(2*pi*f(g)*t+phase);
    for d=1:N/bits_for_symbol
        symbol_E_Inphase_f=sum(Inphase_f_mem((d-1)*D+1 : d*D).^2);        
        Inphase_f_mem((d-1)*D+1 : d*D) = Inphase_f_mem((d-1)*D+1 : d*D)/sqrt(symbol_E_Inphase_f);      
    end
    Inphase_f(g,:)=Inphase_f_mem(1,:);
end
for l=1:14
%%%%%%%%%%%%%%%%%%%%%%%%%Rx signal with demodulation%%%%%%%%%%%
    noise=randn(1,length(Tx))*sig(l)+1j*randn(1,length(Tx))*sig(l);
    Rx=Tx+noise;
    for g=1:length(f)
        Rx_I(g,:)=real(Rx) .*Inphase_f(g,:);
    end
    Rx_bits=zeros(1,N);
    bit_binary_Index=[];
    for k=1:N/bits_for_symbol
        for g=1:length(f)  
            Rx_I_symbol(g,:)=Rx_I(g,(k-1)*D+1 : k*D);
            z(g)=sum(Rx_I_symbol(g,:));
        end
    %%%%%%%%%%%%%%%decision%%%%%%%%%%%%%%%%%%
        [sup,Index]=max(z); 
        binary_Index=de2bi(Index-1, bits_for_symbol, 'left-msb');
        bit_binary_Index=[bit_binary_Index, binary_Index];       
        if Index==(symbols_decimal(k)+1)
            sym_error_count(l)=sym_error_count(l)+0;
        else 
            sym_error_count(l)=sym_error_count(l)+1;
        end
    end
    bit_error_count(l)=sum(bit_binary_Index ~= bits);
    BER(l)=bit_error_count(l)/N;
    SER(l)=sym_error_count(l)/(N/bits_for_symbol);
end

x=-3:10;
[ber,ser]=berawgn(ebn0_db,'fsk',M,'coherent'); 
semilogy(x,BER, x, ber, '-ro',x, SER,x,ser,'-o');
legend('AWGN BER', 'theo AWGN BER', 'AWGN SER', 'theo AWGN SER');



