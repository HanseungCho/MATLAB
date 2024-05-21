clc;
clear;
T=1/(10^5);
t=-10*T:T/10:10*T; %샘플링 10^6
Fs=10^6; 
b=[0, 0.5, 1];
%RC filter
figure(1)
for n=1:3
    h1_t=sinc(t/T).*cos(pi*b(n)*t/T)./(1-(2*b(n)*t/T).^2);  
    h1_f=fftshift(fft(h1_t));

    f1_axis=linspace(-Fs/2,Fs/2,length(h1_f));
    subplot(2,1,1);
    plot(t,h1_t);
    hold on  
    legend('beta=0','beta=0.5','beta=1');
    xlim([-0.00006, 0.00006]);
    subplot(2,1,2);
    plot(f1_axis,abs(h1_f));
    hold on
    legend('beta=0','beta=0.5','beta=1');
    xlim([-Fs/4,Fs/4])
end
%RRC filter
figure(2)
for n=1:3
    h2_t=(sin(pi*t/T*(1-b(n)))+4*b(n)*t/T.*cos(pi*t/T*(1+b(n))))./(pi*t/T.*(1-(4*b(n)*t/T).^2));  
    h2_t(101)=1+b(n)*(4/pi-1);
    h2_f=fftshift(fft(h2_t));

    f2_axis=linspace(-Fs/2,Fs/2,length(h2_f));
    subplot(2,1,1);
    plot(t,h2_t);
    hold on  
    legend('beta=0','beta=0.5','beta=1');
    xlim([-0.00006, 0.00006]);
    subplot(2,1,2);
    plot(f2_axis,abs(h2_f));
    hold on
    legend('beta=0','beta=0.5','beta=1');
    xlim([-Fs/4,Fs/4])
end
%oversampling
figure(3)
data=randi([0,1],1,1000)*2-1;
oversample_data=zeros(1,10000);
for k=1:1000
oversample_data(10*(k-1)+1)=data(k);
    for l=2:10
        oversample_data(10*(k-1)+l)=0;
    end
end
space=linspace(0,0.01,10000); 
%Fs=10^6로 샘플링된 신호의 sample T=1/(10^6) conv 신호개수 10000개에 fitting
g1=conv(oversample_data,h1_t,'same');
g2=conv(oversample_data,h2_t,'same');
plot(space,oversample_data);
hold on
plot(space,g1);
hold on
plot(space,g2);
legend('oversampling','RC','RRC');
xlim([0 0.0003]);


