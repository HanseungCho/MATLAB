clc; clear;
SR=[1 1 1 1 1 1];
%Polynomial=x^6+x+1;

PN(1)=SR(6);
for i=1:2^length(SR)-2
   mem=mod(SR(6)+SR(1),2);
   SR=circshift(SR,1);
   SR(1)=mem;
   PN(i+1)=SR(6);
end
PNcode=2*PN-1;
shift=0:128;
for i=1:length(shift)
    autocorr(i)=sum(PNcode.*circshift(PNcode,shift(i)));
end
figure(1)
stem(shift,autocorr/length(PNcode)) %PNcode 길이로 정규화
title("PNcode의 Autocorrelation")
%run property
l=1;
run=zeros(1,length(SR));
PNcode(length(PNcode)+1)=2;
for i=2:length(PN)+1   
    if PNcode(i) == PNcode(i-1)
        l=l+1;
    else
        run(l)=run(l)+1;
        l=1;
    end
end
PNcode=PNcode(1:length(PNcode)-1);
%balance property
num1=sum(PN==1);
num0=sum(PN==0);
figure(2)
subplot(2,1,1)
stem(0:1,[num0 num1])
title("PNcode의 Balance property")
xlim([-1 2])
subplot(2,1,2)
stem(run/sum(run))
title("PNcode의 Run property")
%shift-and-add property (shift=3)
shifted_PN=circshift(PN,3);
shifted_add_PN=bitxor(PN,shifted_PN);
for i=0:length(PN)
    mem=circshift(shifted_add_PN,i);
    if mem == PN
        same_shift=i;
    end
end
figure(3)
subplot(3,1,1)
stem(PN)
title("Original PN code")
subplot(3,1,2)
stem(shifted_add_PN)
title("Shift-and-add PN code(shift=3)")
subplot(3,1,3)
stem(circshift(shifted_add_PN,same_shift));
title("Shift-and-add PN code(shift=3) with 31 shift")

%Gold code Simulation
pref_SR=[1 1 1 1 1 1];
%Polynomial=x^6+x^5+x^2+x+1;

pref_PN(1)=pref_SR(6);
for i=1:2^length(pref_SR)-2
   mem=mod(pref_SR(6)+pref_SR(5)+pref_SR(2)+pref_SR(1),2);
   pref_SR=circshift(pref_SR,1);
   pref_SR(1)=mem;
   pref_PN(i+1)=pref_SR(6);
end
pref_PNcode=2*pref_PN-1;
for i=1:length(PNcode)+1
    pref_autocorr(i)=sum(PNcode.*circshift(pref_PNcode,i-1));
end
figure(4)
stem(0:length(pref_autocorr)-1,pref_autocorr/length(pref_PN)) %PNcode 길이로 정규화
title("preferred PN code Autocorrelation(x^6+x+1, x^6+x^5+x^2+x+1), Crosscorrelation(-1/63(-0.015), -17/63(-0.269), 15/63(0.238))")

Gold(1,:)=PN;
Gold(2,:)=pref_PN;
for i=1:length(PNcode)
    Gold(i+2,:)=mod(PN+circshift(pref_PN,i-1),2); 
end
Goldcode=2*Gold-1;

for i=1:length(PNcode)
    Gold_autocorr(i)=sum(Goldcode(6,:).*circshift(Goldcode(6,:),i-1));
end
figure(5)
stem(0:length(Gold_autocorr)-1,Gold_autocorr/length(Gold_autocorr)) %PNcode 길이로 정규화
title("Gold code Autocorrelation (-1/63(-0.015), -17/63(-0.269), 15/63(0.238), 63/63(1))")

for i=1:length(PNcode)
    Gold_crosscorr(i)=sum(Goldcode(6,:).*circshift(Goldcode(10,:),i-1));
end
figure(6)
stem(0:length(Gold_crosscorr)-1,Gold_crosscorr/length(Gold_crosscorr)) %PNcode 길이로 정규화
title("Gold code Crosscorrelation (-1/63(-0.015), -17/63(-0.269), 15/63(0.238))")

%Gold code added with 0
for i=1:length(PNcode)
    Gold_added0_autocorr(i)=sum([Goldcode(6,:) 0].*circshift([Goldcode(6,:) 0],i-1));
end
figure(7)
stem(0:length(Gold_added0_autocorr)-1,Gold_added0_autocorr/length(Gold_added0_autocorr)) %PNcode 길이로 정규화
title("Gold code added with 0 Autocorrelation (-1/63(-0.015), -17/63(-0.269), 15/63(0.238), 63/63(1))")

for i=1:length(PNcode)
    Gold_added0_crosscorr(i)=sum([Goldcode(6,:) 0].*circshift([Goldcode(10,:) 0],i-1));
end
figure(8)
stem(0:length(Gold_added0_crosscorr)-1,Gold_added0_crosscorr/length(Gold_added0_crosscorr)) %PNcode 길이로 정규화
title("Gold code added with 0 Crosscorrelation (-1/63(-0.015), -17/63(-0.269), 15/63(0.238))")