clc; clear;
data_length=800;
cp_length=20;
% -1 or 1 generation in frequency domain
data=2*randi([0,1],1,data_length)-1;
data_fdomain=[zeros(1,112),data,zeros(1,112)];
data_tdomain=ifft(data_fdomain);
Tx=zeros(1,1044);
Tx(1:20)=data_tdomain(1005:1024); %CP 20
Tx(21:1044)=data_tdomain; 
Rx=Tx; %through flat channel
signal_tdomain=Rx(21:1044);
signal_fdomain=fft(signal_tdomain);
signal_fdomain=round(real(signal_fdomain(113:912)));
subplot(4,1,1);
plot(data); title("original data");
subplot(4,1,2);
plot(abs(data_tdomain)); title("original data IFFT");
subplot(4,1,3);
plot(abs(Tx)); title("original data IFFT with CP");
subplot(4,1,4);
plot(signal_fdomain); title("decoded data");
disp(append("errors : ", string(sum(data~=signal_fdomain))));
