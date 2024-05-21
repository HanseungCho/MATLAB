clear
N = 10^5; % number of bits or symbols
SD = 10; % symbol duration
t = 0:1/SD:1-1/SD;

for ia = 1:N
    tR(SD*(ia-1)+1:SD*ia) = t;
end

f1 = 3000;
f2 = 5050;

bitstream = randi([0,1], 1, N);

for n = 1:N
    if bitstream(n) == 0
        freqdum(n) = f1;
    elseif bitstream(n) == 1
        freqdum(n) = f2;
    end
end


for ib = 1:N
    freq(SD*(ib-1)+1:SD*ib) = freqdum(ib);
end

x = (sqrt(2/SD))*cos(2*pi*freq.*tR); %FSK signal 
n = (randn(1,N*SD)+rand(1,N*SD)*j)*sqrt(1/2); % white gaussian noise

ebn0 = -3:10;
EN=10.^(ebn0/10);

for j = 1:14
    k = 0.5;
    sig(j) = (1/(2*k*EN(j)))^(1/2);
end

for ic = 1:14
    gn = n*sig(ic);
    y = x + gn;

    y1 = real(y.*(sqrt(2/SD)*cos(2*pi*f1.*tR)));
    y2 = real(y.*(sqrt(2/SD)*cos(2*pi*f2.*tR)));

    for fac = 1:N
        itg1(fac) = sum(y1((fac-1)*SD+1:fac*SD));
        itg2(fac) = sum(y2((fac-1)*SD+1:fac*SD));
    end

    comp12 = itg1(1:1:end) < itg2(1:1:end);

    a = 0;
    for u = 1:N
        if bitstream(u) ~= comp12(u)
            a = a+1;
        end
    end
    
    BER(ic) = a/N;
end

ebn0_db=-3:10; 
ebn0=zeros(1,14);
for k=1:14
ebn0(k)=10^(ebn0_db(k)/10);
sig(k)=sqrt(1/(2*ebn0(k)));
end
ber = 0.5*erfc(sqrt((10.^(ebn0/10))/2)); %theoretical BER 
ser = ber/N;
x=-3:10;
[ber0,ser0]=berawgn(ebn0_db,'fsk',2,'coherent'); 
figure
semilogy(ebn0_db,ber0);
hold on
semilogy(ebn0_db,BER);
axis([-3 10 10^-4 0.5])
grid on
legend('theoBER', 'BER');
xlabel('Eb/No(dB)')
ylabel('Error rate')
title('BER')