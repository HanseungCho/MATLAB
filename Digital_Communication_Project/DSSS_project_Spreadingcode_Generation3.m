%Tc=1/(10^7);
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
%딜레이 추가
delay=randi([1,length(JPL_C)]);
Delayed_JPL_C=circshift(JPL_C, delay);
Delayed_K5_C=circshift(K5_C, delay);

