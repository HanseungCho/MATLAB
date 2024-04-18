EsN0=1:1:10; %10dB
M=8; 
syms k M
SER1_func=exp(-EsN0)*symsum(((-1)^(k+1))*nchoosek(M-1,k+1)*1/(M-(k+1))*(exp(EsN0)).^(1/(k+1)),k,1,M-1);
