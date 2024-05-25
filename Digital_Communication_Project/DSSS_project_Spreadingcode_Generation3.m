%Tc=1/(10^7);
%JPL & K5 generation
clc; clear; close all
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
%Local Generated Code
rep_C1=repmat(C1, 1, length(JPL_C)/length(C1));
rep_C2=repmat(C2, 1, length(JPL_C)/length(C2));
rep_C3=repmat(C3, 1, length(JPL_C)/length(C3));
rep_C4=repmat(C4, 1, length(JPL_C)/length(C4));
rep_C5=repmat(C5, 1, length(JPL_C)/length(C5));
snr=-15:10;
%딜레이 추가

N=1000;
for m=1:length(snr)
    error=0;
    for n=1:N
        delay=randi([1,length(JPL_C)]);
        Delayed_JPL_C=circshift(JPL_C, delay);
        Delayed_K5_C=circshift(K5_C, delay);
        Original_JPL_C=Delayed_JPL_C;
        Original_K5_C=Delayed_K5_C;
        Delayed_JPL_C=awgn(Delayed_JPL_C, snr(m));
        Delayed_K5_C=awgn(Delayed_K5_C, snr(m));
        
        %Autocorrelation
        shift=0:43890*2+1;
        for i=1:length(shift)
            JPLC1_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C1,i-1));
        end
        for i=1:length(shift)
            JPLC2_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C2,i-1));
        end
        for i=1:length(shift)
            JPLC3_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C3,i-1));
        end
        for i=1:length(shift)
            JPLC4_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C4,i-1));
        end
        for i=1:length(shift)
            JPLC5_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C5,i-1));
        end
        [v1,index1]=max(JPLC1_autocorr(1:length(C1)));
        [v2,index2]=max(JPLC2_autocorr(1:length(C2)));
        [v3,index3]=max(JPLC3_autocorr(1:length(C3)));
        [v4,index4]=max(JPLC4_autocorr(1:length(C4)));
        [v5,index5]=max(JPLC5_autocorr(1:length(C5)));
        index1=index1-1;%index에서 1을 빼면 각 코드의 timing offset
        index2=index2-1;
        index3=index3-1;
        index4=index4-1;
        index5=index5-1;
        Offset_D1=circshift(D1,index1);
        Offset_D2=circshift(D2,index2);
        Offset_D3=circshift(D3,index3);
        Offset_D4=circshift(D4,index4);
        Offset_D5=circshift(D5,index5);
        for i=1:43890
            if (Offset_D1(1)+Offset_D2(1)+Offset_D3(1)+Offset_D4(1)+Offset_D5(1)) >= 3
                Regen_JPL_D(i)=1;
            else
                Regen_JPL_D(i)=0;
            end
            Offset_D1=circshift(Offset_D1,-1);
            Offset_D2=circshift(Offset_D2,-1);
            Offset_D3=circshift(Offset_D3,-1);
            Offset_D4=circshift(Offset_D4,-1);
            Offset_D5=circshift(Offset_D5,-1);
        end
        Regen_JPL_C=2*Regen_JPL_D-1;
        
        if sum(Original_JPL_C ~= Regen_JPL_C) == 0
    
        else 
            error=error+1;
        end
    end
    SP(m)=error/N;
end
plot(snr,SP)

