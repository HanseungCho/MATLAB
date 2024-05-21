clc;
clear;
o=randi([0,1],1,1000);
x0=zeros(1,100000);
x1=zeros(1,100000);
x2=zeros(1,100000);
x3=zeros(1,100000);
space=linspace(0,1,100000);

for i0 =1:1000
    x0(100*(i0-1)+1)=o(i0);
    for c0 = 2:100
        k0=100*(i0-1)+c0;
        x0(k0)=0;
    end
end

for i1 = 1:1000
    for c1 = 1:100
        k1=100*(i1-1)+c1;
        x1(k1)=(o(i1)-1/2)*2;
    end
end

for i2 = 1:1000
    for c2 = 1:50
        if o(i2)==0
        k2=100*(i2-1)+c2;
        x2(k2)=-1;
        x2(k2+50)=0;
        else
        k2=100*(i2-1)+c2;
        x2(k2)=1;
        x2(k2+50)=0;
        end
    end
end

for i3 = 1:1000
    for c3 = 1:50
        if o(i3)==0
        k3=100*(i3-1)+c3;
        x3(k3)=-1;
        x3(k3+50)=1;
        else
        k3=100*(i3-1)+c3;
        x3(k3)=1;
        x3(k3+50)=-1;
        end
    end
end
figure(1)
subplot(4,1,1);
stem([1:10],o(1:10))
subplot(4,1,2);
plot(space(1:1000),x1(1:1000))
subplot(4,1,3);
plot(space(1:1000),x2(1:1000))
subplot(4,1,4);
plot(space(1:1000),x3(1:1000))

figure(2)
L=100000;
fs=100000;
t = 0:(1/fs):(10-1/fs);
f=linspace(-fs/2,fs/2,L);
y1=fft(x1);
y11=abs(fftshift(y1));
y2=fft(x2);
y22=abs(fftshift(y2));
y3=fft(x3);
y33=abs(fftshift(y3));
plot(f,y11)
hold on
plot(f,y22)
hold on
plot(f,y33)
xlim([-2000 2000])
ylim([0 9000])
