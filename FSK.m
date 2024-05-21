function [symbols ,Rx0, Tx, noise]=FSK(M,N,T_bit,samplingrate_bit,Ebn0_db,coherence)
ebn0_db=-3:10; 
ebn0=zeros(1,14);
sig=zeros(1,14);
phase=pi/7;
for k=1:14
ebn0(k)=10^(ebn0_db(k)/10); 
sig(k)=sqrt(1/(2*ebn0(k)));
end
bits=randi([0,1],1,N);
symbols=zeros(1,2);
if (M==4)
    for l=1:N/(M/2)
        symbols(l)=2*bits(2*l-1)+bits(2*l);
    end
else
    for l=1:N/(M/2)
        symbols(l)=bits(l);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (M == 4) && (coherence == 0) % 4FSK noncoherent
    T=2*T_bit;
    D=2*samplingrate_bit; %D samplingrte_symbol
    Fs=1/T*D;
    t=0:(1/Fs):D*(N/sqrt(M))*(1/Fs)-(1/Fs);
    f1=Fs/50;
    f2=Fs/50+1/(T);
    f3=Fs/50+2/(T);
    f4=Fs/50+3/(T);
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
    noise=randn(1,length(Tx))*sig(Ebn0_db+4);
    Rx0=Tx+noise;
elseif (M == 4) && (coherence == 1) % 4FSK coherent
    T=2*T_bit;
    D=2*samplingrate_bit;
    Fs=1/T*D;
    t=0:(1/Fs):D*(N/sqrt(M))*(1/Fs)-(1/Fs);
    t_symbol=0:(1/Fs):D*(1/Fs)-(1/Fs);
    coherent_t=zeros(1,N/sqrt(M)*D);
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
    noise=randn(1,length(Tx))*sig(Ebn0_db+4);
    Rx0=Tx+noise;

elseif (M == 2) && (coherence == 0)
    T=T_bit;
    D=samplingrate_bit;
    Fs=1/T*D;
    t=0:(1/Fs):D*N*(1/Fs)-(1/Fs);
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
    noise=randn(1,length(Tx))*sig(Ebn0_db+4);
    Rx0=Tx+noise;
else
    T=T_bit;
    D=samplingrate_bit;
    Fs=1/T*D;
    bits=randi([0,1],1,N);
    t=0:(1/Fs):D*N*(1/Fs)-(1/Fs);
    t_symbol=0:(1/Fs):D*(1/Fs)-(1/Fs);
    coherent_t=zeros(1,N*D);
    for v=1:N
        coherent_t((v-1)* D +1 : v * D) = t(1 : D);
    end
    f1=Fs/50;
    f2=Fs/50+1/(2*T);
    symbol_f1=cos(2*pi*f1*t_symbol+phase);
    E_f1=sum(symbol_f1.^2);
    symbol_f1=symbol_f1/sqrt(E_f1); %symbol energy를 1로 맞춤.
    symbol_f2=cos(2*pi*f2*t_symbol+phase);
    E_f2=sum(symbol_f2.^2);
    symbol_f2=symbol_f2/sqrt(E_f2);
    Tx = [];
    for b=1:N
        if bits(b) == 0
            Tx = [Tx symbol_f1];
        else
            Tx = [Tx symbol_f2];
        end
    end
    noise=randn(1,length(Tx))*sig(Ebn0_db+4);
    Rx0=Tx+noise;
end
end
