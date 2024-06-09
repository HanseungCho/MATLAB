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
    K5_C(i)=sign(3*C1(1)+C2(1)-C3(1)-C4(1)+C5(1));
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
snr=-35:-12;

N=10^4;
for m=1:length(snr) %SNR에 따라 성능 검증
    JPL_detection=0;
    K5_detection=0;
    for n=1:N
        %랜덤한 딜레이를 추가하여 신호 생성
        delay=randi([1,length(JPL_C)]);
        Delayed_JPL_C=circshift(JPL_C, delay);
        Delayed_K5_C=circshift(K5_C, delay);
        Original_JPL_C=Delayed_JPL_C;
        Original_K5_C=Delayed_K5_C;
        %신호에 SNR에 따른 AWGN 추가
        Delayed_JPL_C=awgn(Delayed_JPL_C, snr(m));
        Delayed_K5_C=awgn(Delayed_K5_C, snr(m));
        
        %각 branch 별로 Autocorrelation
        for i=1:length(C1)
            JPLC1_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C1,i-1));
        end
        for i=1:length(C2)
            JPLC2_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C2,i-1));
        end
        for i=1:length(C3)
            JPLC3_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C3,i-1));
        end
        for i=1:length(C4)
            JPLC4_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C4,i-1));
        end
        for i=1:length(C5)
            JPLC5_autocorr(i)=sum(Delayed_JPL_C.*circshift(rep_C5,i-1));
        end
        %각 branch 별로 autocorrelation 피크 지점 탐지로 phase 확보
        [v1,index1]=max(JPLC1_autocorr(1:length(C1)));
        [v2,index2]=max(JPLC2_autocorr(1:length(C2)));
        [v3,index3]=max(JPLC3_autocorr(1:length(C3)));
        [v4,index4]=max(JPLC4_autocorr(1:length(C4)));
        [v5,index5]=max(JPLC5_autocorr(1:length(C5)));
        index1=index1-1;%index에서 1을 빼서 0점 보정
        index2=index2-1;
        index3=index3-1;
        index4=index4-1;
        index5=index5-1;
        %앞에서 찾은 phase만큼 original C1~C5 코드를 딜레이 시킨후 JPL코드 합성
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
            JPL_detection=JPL_detection+1;  
        else 
            
        end
        %K5 Autocorrelation
        %각 branch 별로 Autocorrelation
        for i=1:length(C1)
            K5C1_autocorr(i)=sum(Delayed_K5_C.*circshift(rep_C1,i-1));
        end
        for i=1:length(C2)
            K5C2_autocorr(i)=sum(Delayed_K5_C.*circshift(rep_C2,i-1));
        end
        for i=1:length(C3)
            K5C3_autocorr(i)=sum(Delayed_K5_C.*circshift(rep_C3,i-1));
        end
        for i=1:length(C4)
            K5C4_autocorr(i)=sum(Delayed_K5_C.*circshift(rep_C4,i-1));
        end
        for i=1:length(C5)
            K5C5_autocorr(i)=sum(Delayed_K5_C.*circshift(rep_C5,i-1));
        end
        %각 branch 별로 autocorrelation 피크 지점 탐지로 phase 확보
        [K5v1,K5index1]=max(K5C1_autocorr(1:length(C1)));
        [K5v2,K5index2]=max(K5C2_autocorr(1:length(C2)));
        [K5v3,K5index3]=max(abs(K5C3_autocorr(1:length(C3))));
        [K5v4,K5index4]=max(abs(K5C4_autocorr(1:length(C4))));
        [K5v5,K5index5]=max(K5C5_autocorr(1:length(C5)));
        K5index1=K5index1-1;%index에서 1을 빼면 각 코드의 timing offset
        K5index2=K5index2-1;
        K5index3=K5index3-1;
        K5index4=K5index4-1;
        K5index5=K5index5-1;
        K5Offset_C1=circshift(C1,K5index1);
        K5Offset_C2=circshift(C2,K5index2);
        K5Offset_C3=circshift(C3,K5index3);
        K5Offset_C4=circshift(C4,K5index4);
        K5Offset_C5=circshift(C5,K5index5);
        %앞에서 찾은 phase만큼 original C1~C5 코드를 딜레이 시킨후 K5코드 합성
        for i=1:43890
            Regen_K5_C(i)=sign(3*K5Offset_C1(1)+K5Offset_C2(1)-K5Offset_C3(1)-K5Offset_C4(1)+K5Offset_C5(1));
            K5Offset_C1=circshift(K5Offset_C1,-1);
            K5Offset_C2=circshift(K5Offset_C2,-1);
            K5Offset_C3=circshift(K5Offset_C3,-1);
            K5Offset_C4=circshift(K5Offset_C4,-1);
            K5Offset_C5=circshift(K5Offset_C5,-1);
        end
        if sum(Original_K5_C ~= Regen_K5_C) == 0
            K5_detection=K5_detection+1;
        else 
            
        end
        
    end
    JPL_SP(m)=JPL_detection/N;
    K5_SP(m)=K5_detection/N;
end
semilogy(snr,JPL_SP)
hold on
semilogy(snr,K5_SP)
legend("JPL Detection probability", "K5 Detection probability")