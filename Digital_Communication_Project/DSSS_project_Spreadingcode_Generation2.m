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