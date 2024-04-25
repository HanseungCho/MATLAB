N_bit=10^2;
bits=2*randi([0,1], N_bit,1)-1; %BPSK
%임의의 PSK 신호 생성
%BPSK
M1=2;
symbol1 = randi([0 M1-1], N_bit, 1);
s1=pskmod(symbol1, M1);
%QPSK
M2=4;
symbol2 = randi([0 M2-1], N_bit, 1);
s2=pskmod(symbol2, M2, pi/M1);
oversamplingrate=100; %1chip=100samples
PN=comm.PNSequence('Polynomial',[1 0 0 0 1 1 0 1], 'SamplesPerFrame', 100, 'InitialConditions',[0 0 0 0 0 0 1]);
pn=PN();
Processing_Gain=length(pn);
oversampled_s1_bits=repelem(s1, oversamplingrate*Processing_Gain);
oversampled_s2_bits=repelem(s2, oversamplingrate*Processing_Gain);
ebn0_db=-3:10; 
esn0=zeros(1,14);
sig=zeros(1,14);
bits_for_symbol=1;
for k=1:14
    esn0_db=ebn0_db+10*log10(bits_for_symbol); 
    esn0(k)=10^(esn0_db(k)/10); 
    sig(k)=sqrt(1/(2*esn0(k)));
end

t=linspace(0,(length(bits)*Tb)-(1/fs),length(oversampled_spreaded_bits));
f=linspace(-fs/2,fs/2,length(oversampled_spreaded_bits));

mem1=cos(2*pi*fc*0.8*t+angle(oversampled_s1_bits.'));
mem2=cos(2*pi*fc*1.05*t+angle(oversampled_s2_bits.'));
Tx_s1=[];
    Tx_s2=[];
    for p=1:length(bits)*length(pn)
        mem_s1=mem1(oversamplingrate*(p-1)+1:oversamplingrate*p);
        E=sum(mem_s1.^2);
        Tx_s1=[Tx_s1 mem_s1/sqrt(E)];
        mem_s2=mem2(oversamplingrate*(p-1)+1:oversamplingrate*p);
        E=sum(mem_s2.^2);
        Tx_s2=[Tx_s2 mem_s2/sqrt(E)];
    end    
for l=1:14    
    noise=randn(1,length(Tx))*sig(l)+1j*randn(1,length(Tx))*sig(l);
    Rx=Tx+noise;

