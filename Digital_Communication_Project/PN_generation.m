function [PN, PNcode]=PN_generation(n, polynomial)
    SR=ones(1,n);
    PN(1)=SR(n);
    for i=1:2^length(SR)-2
       mem=mod(sum(SR(polynomial)),2);
       SR=circshift(SR,1);
       SR(1)=mem;
       PN(i+1)=SR(n);
    end
    PNcode=2*PN-1;
end