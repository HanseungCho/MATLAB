fs=10000;
t=linspace(0,1,10000);
f=linspace(-fs/2,fs/2,10000);
x=randn(1,10000);

figure(1)
subplot(4,1,1)
plot(t,x)
title("noise")
subplot(4,1,2)
plot(f, abs(fftshift(fft(x))))
title("noise fft")
subplot(4,1,3)
plot(t, x.^2);
title("noise^2")
subplot(4,1,4)
plot(f, abs(fftshift(fft(x.^2))));
title("noise^2 fft")

x_LPF=lowpass(x, 500, fs);
figure(2)
subplot(4,1,1)
plot(t, x_LPF)
title("LPF filterd noise")
subplot(4,1,2)
plot(f, abs(fftshift(fft(x_LPF))))
title("LPF filterd noise fft")
subplot(4,1,3)
t_corr=linspace(-1, 1, length(xcorr(x)));
plot(t_corr, xcorr(x))
title("noise autocorrelation")
subplot(4,1,4)
plot(t_corr, xcorr(x_LPF))
title("LPF filtered noise autocorrelation")
ylim([0 150])

figure(3)
bpsk=2*randi([0 1],1,1000)-1;
oversampling_rate=10;
for k=1:length(bpsk)
    PAM(oversampling_rate*(k-1)+1:oversampling_rate*k)=bpsk(k);
end
subplot(3,1,1)
plot(t, PAM);
title("PAM")
xlim([0 0.1])
subplot(3,1,2)
plot(t_corr, xcorr(PAM));
title("autocorrelation of PAM")
subplot(3,1,3)

tau=(0:1:(10000-1))/10000;
taau=uint32(tau*10000);
PAM_quadratic=PAM.*(circshift(PAM, taau));
cyclic_autocorr=abs(fftshift(fft(PAM_quadratic)));
cyclic_f=linspace(-fs/2,fs/2,length(cyclic_autocorr));
plot3(cyclic_f, tau, cyclic_autocorr);
c
title('PAM cyclic autocorrelation')
