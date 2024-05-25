clc; clear; close all
%JPL & K5 generation
C1=[1 -1];
C2=[1 1 1 -1 -1 1 -1];
C3=[1 1 1 -1 -1 -1 1 -1 1 1 -1];
C4=[1 1 1 1 -1 -1 -1 1 -1 -1 1 1 -1 1 -1];
C5=[1 1 1 1 -1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 -1];
D1=(C1+1)/2;
D2=(C2+1)/2;
D3=(C3+1)/2;
D4=(C4+1)/2;
D5=(C5+1)/2;

for i=1:43890
    if (D1(1)+D2(1)+D3(1)+D4(1)+D5(1)) >= 3
        JPL_D(i)=1;
    else
        JPL_D(i)=0;
    end
    D1=circshift(D1,-1);
    D2=circshift(D2,-1);
    D3=circshift(D3,-1);
    D4=circshift(D4,-1);
    D5=circshift(D5,-1);
end
JPL_C=2*JPL_D-1;
for i=1:43890
    K5_C(i)=sign(3*C1(1)-C2(1)+C3(1)+C4(1)-C5(1));
    C1=circshift(C1,-1);
    C2=circshift(C2,-1);
    C3=circshift(C3,-1);
    C4=circshift(C4,-1);
    C5=circshift(C5,-1);
end
K5_D=(K5_C+1)/2;
%Autocorrelation
shift=0:43890*2+1;
for i=1:length(shift)
    JPL_autocorr(i)=sum(JPL_C.*circshift(JPL_C,shift(i)));
end
for i=1:length(shift)
    K5_autocorr(i)=sum(K5_C.*circshift(K5_C,shift(i)));
end
figure(1)
subplot(2,1,1)
stem(shift,JPL_autocorr/length(JPL_C)) %PNcode 길이로 정규화
title("JPL Autocorrelation")
subplot(2,1,2)
stem(shift,K5_autocorr/length(K5_C)) %PNcode 길이로 정규화
title("K5 Autocorrelation")

%transition per chip
l=0;
for i=2:length(JPL_D)   
    if JPL_D(i) ~= JPL_D(i-1)
        l=l+1;
    else
    end
end
JPL_TPC=l/length(JPL_D);
p=0;
for i=2:length(K5_D)   
    if K5_D(i) ~= K5_D(i-1)
        p=p+1;
    else
    end
end
K5_TPC=p/length(K5_D);
%DC bias
JPL_DC=mean(JPL_C);
K5_DC=mean(K5_C);

%IS-95 CDMA walsh code generation
N=15; %2^15=32,768 length walsh code
H=1;
for i=1:N
    H=[H H; H -H];
end
Walsh_C=H(2,:);
Walsh_D=(Walsh_C+1)/2;

shift=0:length(Walsh_D)*2+1;
for i=1:length(shift)
    Walsh_autocorr(i)=sum(Walsh_C.*circshift(Walsh_C,shift(i)));
end
figure(2)
stem(shift,Walsh_autocorr/length(Walsh_C)) %PNcode 길이로 정규화
title("IS-95 Walsh Autocorrelation")
%transition per chip
k=0;
for i=2:length(Walsh_D)   
    if Walsh_D(i) ~= Walsh_D(i-1)
        k=k+1;
    else
    end
end
Walsh_TPC=k/length(Walsh_D);
%DC bias
Walsh_DC=mean(Walsh_C);