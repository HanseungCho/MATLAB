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

%IS-95 PN code
%Inphase PN code polynomial: x^15+x^13+x^9+x^8+x^7+x^5+1
SRI=ones(1,15);
PNI(1)=SRI(15);
for i=1:2^length(SRI)-2
   memI=mod(SRI(15)+SRI(13)+SRI(9)+SRI(8)+SRI(7)+SRI(5),2);
   SRI=circshift(SRI,1);
   SRI(1)=memI;
   PNI(i+1)=SRI(6);
end
PNcodeI=2*PNI-1;
%Quadrature PN code polynomial: x^15+x^12+x^11+x^10+x^6+x^5+x^4+x^3+1
SRQ=ones(1,15);
PNQ(1)=SRQ(15);
for i=1:2^length(SRQ)-2
   memQ=mod(SRQ(15)+SRQ(12)+SRQ(11)+SRQ(10)+SRQ(6)+SRQ(5)+SRQ(4)+SRQ(3),2);
   SRQ=circshift(SRQ,1);
   SRQ(1)=memQ;
   PNQ(i+1)=SRQ(6);
end
PNcodeQ=2*PNQ-1;
%Autocorrelation
shift=0:length(PNcodeI)+1;
for i=1:length(shift)
    PNcodeI_autocorr(i)=sum(PNcodeI.*circshift(PNcodeI,shift(i)));
end
figure(2) 
plot(shift, PNcodeI_autocorr/length(PNcodeI_autocorr))
title("PNcode(Inphase) Autocorrelation")
shift=0:length(PNcodeQ)+1;
for i=1:length(shift)
    PNcodeQ_autocorr(i)=sum(PNcodeQ.*circshift(PNcodeQ,shift(i)));
end

figure(3) 
plot(shift, PNcodeQ_autocorr/length(PNcodeQ_autocorr))
title("PNcode(Quadrature) Autocorrelation")
%transition per chip
l=0;
for i=2:length(PNcodeI)   
    if PNcodeI(i) ~= PNcodeI(i-1)
        l=l+1;
    else
    end
end
PNcodeI_TPC=l/length(PNcodeI);
p=0;
for i=2:length(PNcodeQ)   
    if PNcodeQ(i) ~= PNcodeQ(i-1)
        p=p+1;
    else
    end
end
PNcodeQ_TPC=p/length(PNcodeQ);
%DC bias
PNcodeI_DC=mean(PNcodeI);
PNcodeQ_DC=mean(PNcodeQ);
