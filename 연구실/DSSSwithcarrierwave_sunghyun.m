clear;
clc;

N = 100000; % 심볼의 갯수
M = 16; % M-ary
k = log2(M); % 심볼당 비트 수

T = 10^(-3); % 심볼듀레이션
samplepersymbol = 100; % 심볼당 샘플링 수
fractiont = [0:T/samplepersymbol:T-T/samplepersymbol];
wholet = [0:T/samplepersymbol:N*T-T/samplepersymbol];
fc = 50000;

%%%% spreading code
spreadingsequence = randi([0,1],1,samplepersymbol);
spreadingcode= 2*spreadingsequence -1; % 1과 -1로 치환
seqlength = length(spreadingcode);

%%%% 16 QAM modulation
Re = [-(2*sqrt(M)/2-1):2:2*sqrt(M)/2-1];
Im = [-(2*sqrt(M)/2-1):2:2*sqrt(M)/2-1];

EbN0_dB = [0:1:15];
EsN0_dB = EbN0_dB + 10*log10(k); % k*EbN0 = EsN0

%%%% Gray code
map = [0,1,3,2];

bitstream = randi([0,1],1,N*k);
bitsmatrix = zeros(N,4);
for i = 1:N
    bitsmatrix(i,1) = bitstream(k*i-3);
    bitsmatrix(i,2) = bitstream(k*i-2);
    bitsmatrix(i,3) = bitstream(k*i-1);
    bitsmatrix(i,4) = bitstream(k*i);
end

bin2DecMatrix = ones(N,1)*(2.^[(k/2-1):-1:0]);

%%%% real
reBitMatrix = bitsmatrix(:,[1:2:k]);
reDecMatrix = sum(reBitMatrix.*bin2DecMatrix,2);

%%%% imaginary
imagBitMatrix = bitsmatrix(:,[2:2:k]);
imagDecMatrix = sum(imagBitMatrix.*bin2DecMatrix,2);

%%%% mapping Gray coded symbols into constellation
modRe = zeros(N,1);
for a = 1:N
    if reDecMatrix(a) == 0
        modRe(a) = -3;
    elseif reDecMatrix(a) == 1
        modRe(a) = -1;
    elseif reDecMatrix(a) == 2
        modRe(a) = 3;
    elseif reDecMatrix(a) == 3
        modRe(a) = 1;
    end
end

modIm = zeros(N,1);
for a = 1:N
    if imagDecMatrix(a) == 0
        modIm(a) = -3;
    elseif imagDecMatrix(a) == 1
        modIm(a) = -1;
    elseif imagDecMatrix(a) == 2
        modIm(a) = 3;
    elseif imagDecMatrix(a) == 3
        modIm(a) = 1;
    end
end

mod = modRe + j*modIm;
k_16QAM = 1/sqrt(10);

%%%% modulation on carrier wave
pastemod = zeros(N,samplepersymbol);
for i = 1:N
    pastemod(i,:)=mod(i);
end

carrierwavematrix = zeros(N,samplepersymbol);
carrierwave = cos(2*pi*fc*fractiont);
sincarrierwave = sin(2*pi*fc*fractiont);
for i = 1:N
    carrierwavematrix(i,:) = carrierwave;
end

modulated_carrierwave_matrix = k_16QAM*pastemod.*carrierwavematrix; % normalization of transmit power to one

spreadingcodematrix = zeros(N,samplepersymbol);
spreadingcodevector = zeros(1,N*samplepersymbol);
for i = 1:N
    spreadingcodematrix(i,:) = spreadingcode;
    spreadingcodevector(1,samplepersymbol*(i-1)+1:samplepersymbol*i) = spreadingcode;
end

%%%% tx
Tx = modulated_carrierwave_matrix.*spreadingcodematrix;

tx_mod = zeros(1,samplepersymbol*N);
for b = 1:N
    tx_mod(seqlength*(b-1)+1:seqlength*b) = Tx(b,:);
end

%%%% noise
EsN0=10.^(EsN0_dB/10); 
sig=sqrt(seqlength./(2*EsN0));
noise_matrix = sig'*[randn(1,length(tx_mod)) + j*randn(1,length(tx_mod))];

%% without noise
        refsymbol = zeros(M,k);
        for refi = 1:M/2
            refsymbol(2*refi,4) = 1;
        end
        for refi = 1:M/4
            refsymbol(4*(refi-1/2)+1:4*refi,3) = 1;
        end
        for refi = 1:M/8
            refsymbol(8*(refi-1/2)+1:8*refi,2) = 1;
        end
        for refi = 1:M/16
            refsymbol(16*(refi-1/2)+1:16*refi,1) = 1;
        end

        refbin2DecMatrix = ones(M,1)*(2.^[(k/2-1):-1:0]);

        %%%% real
        refreBitMatrix = refsymbol(:,[1:2:k]);
        refreDecMatrix = sum(refreBitMatrix.*refbin2DecMatrix,2);

        %%%% imaginary
        refimagBitMatrix = refsymbol(:,[2:2:k]);
        refimagDecMatrix = sum(refimagBitMatrix.*refbin2DecMatrix,2);

        %%%% mapping Gray coded symbols into constellation
        refmodRe = zeros(M,1);
        for a = 1:M
            if refreDecMatrix(a) == 0
                refmodRe(a) = -3;
            elseif refreDecMatrix(a) == 1
                refmodRe(a) = -1;
            elseif refreDecMatrix(a) == 2
                refmodRe(a) = 3;
            elseif refreDecMatrix(a) == 3
                refmodRe(a) = 1;
            end
        end

        refmodIm = zeros(M,1);
        for a = 1:M
            if refimagDecMatrix(a) == 0
                refmodIm(a) = -3;
            elseif refimagDecMatrix(a) == 1
                refmodIm(a) = -1;
            elseif refimagDecMatrix(a) == 2
                refmodIm(a) = 3;
            elseif refimagDecMatrix(a) == 3
                refmodIm(a) = 1;
            end
        end

        refmod = refmodRe + j*refmodIm;

        %%%% modulation on carrier wave
        refpastemod = zeros(M,samplepersymbol);
        for i = 1:M
            refpastemod(i,:)=refmod(i);
        end

        refcarrierwavematrix = zeros(M,samplepersymbol);
        for i = 1:M
            refcarrierwavematrix(i,:) = carrierwave;
        end

        refmodulated_carrierwave_matrix = k_16QAM*refpastemod.*refcarrierwavematrix; % normalization of transmit power to one

        refspreadingcodematrix = zeros(M,samplepersymbol);
        refspreadingcodevector = zeros(1,M*samplepersymbol);
        for i = 1:M
            refspreadingcodematrix(i,:) = spreadingcode;
            refspreadingcodevector(1,samplepersymbol*(i-1)+1:samplepersymbol*i) = spreadingcode;
        end

        %%%% tx
        ref_tx = refmodulated_carrierwave_matrix.*refspreadingcodematrix;

        ref_tx_mod = zeros(1,samplepersymbol*M);
        for b = 1:M
            ref_tx_mod(seqlength*(b-1)+1:seqlength*b) = ref_tx(b,:);
        end

        %%%% despreading
        ref_spread_ref_tx_mod=refspreadingcodevector.*ref_tx_mod;

        ref_spread_ref_tx_mod_matrix = zeros(M,samplepersymbol);
        for i = 1:M
            ref_spread_ref_tx_mod_matrix(i,:) = ref_spread_ref_tx_mod(1,samplepersymbol*(i-1)+1:samplepersymbol*i);
        end

        %%%% correlator with cosine carrier wave
        ref_outputcorrelator = ref_spread_ref_tx_mod_matrix*(carrierwave)';

        % scatter(real(ref_outputcorrelator),imag(ref_outputcorrelator))

        %%%% criterion
        ax = sort(real(ref_outputcorrelator));
        axcenter = zeros(1,k-1);
        for axi = 1:k-1
            axcenter(axi) = (ax(4*axi)+ax(4*(axi+1)))/2;
        end
%% with noise
%%%% despreading
for ebn0index = 1:length(EbN0_dB)
    noisy_tx_mod = tx_mod + noise_matrix(ebn0index,:);

    %%%% despreading
    spread_noisy_tx_mod=spreadingcodevector.*noisy_tx_mod;

    spread_noisy_tx_mod_matrix = zeros(N,samplepersymbol);
    for i = 1:N
        spread_noisy_tx_mod_matrix(i,:) = spread_noisy_tx_mod(1,samplepersymbol*(i-1)+1:samplepersymbol*i);
    end

    %%%% correlator with cosine carrier wave
    outputcorrelator = spread_noisy_tx_mod_matrix*(carrierwave)';

    %%%% demodulation
    y_Re = real(outputcorrelator);
    y_Im = imag(outputcorrelator);

    demod_Rx_symbol_Re = zeros(1,N);

    for c = 1:N
        if (y_Re(c) >= axcenter(3))
            demod_Rx_symbol_Re(c) = 2;
        elseif (y_Re(c) < axcenter(3))&&(y_Re(c) >= axcenter(2))
            demod_Rx_symbol_Re(c) = 3;
        elseif (y_Re(c) < axcenter(2))&&(y_Re(c) >= axcenter(1))
            demod_Rx_symbol_Re(c) = 1;
        elseif (y_Re(c) < axcenter(1))
            demod_Rx_symbol_Re(c) = 0;
        end
    end

    demod_Rx_symbol_Im = zeros(1,N);

    for c = 1:N
        if (y_Im(c) >= axcenter(3))
            demod_Rx_symbol_Im(c) = 2;
        elseif (y_Im(c) < axcenter(3))&&(y_Im(c) >= axcenter(2))
            demod_Rx_symbol_Im(c) = 3;
        elseif (y_Im(c) < axcenter(2))&&(y_Im(c) >= axcenter(1))
            demod_Rx_symbol_Im(c) = 1;
        elseif (y_Im(c) < axcenter(1))
            demod_Rx_symbol_Im(c) = 0;
        end
    end

    symerror(ebn0index) = 0;
    for i = 1:N
        if (demod_Rx_symbol_Re(i) ~= reDecMatrix(i))|(demod_Rx_symbol_Im(i) ~= imagDecMatrix(i))
            symerror(ebn0index) = symerror(ebn0index) + 1;
        end
    end

    ser = symerror/N;

    demod_Rx_bit = zeros(1,N*k);

    for c = 1:N
        if demod_Rx_symbol_Re(c) == 0
            demod_Rx_bit(4*c-3) = 0;
            demod_Rx_bit(4*c-1) = 0;
        elseif demod_Rx_symbol_Re(c) == 1
            demod_Rx_bit(4*c-3) = 0;
            demod_Rx_bit(4*c-1) = 1;
        elseif demod_Rx_symbol_Re(c) == 2
            demod_Rx_bit(4*c-3) = 1;
            demod_Rx_bit(4*c-1) = 0;
        elseif demod_Rx_symbol_Re(c) == 3
            demod_Rx_bit(4*c-3) = 1;
            demod_Rx_bit(4*c-1) = 1;
        end
        if demod_Rx_symbol_Im(c) == 0
            demod_Rx_bit(4*c-2) = 0;
            demod_Rx_bit(4*c) = 0;
        elseif demod_Rx_symbol_Im(c) == 1
            demod_Rx_bit(4*c-2) = 0;
            demod_Rx_bit(4*c) = 1;
        elseif demod_Rx_symbol_Im(c) == 2
            demod_Rx_bit(4*c-2) = 1;
            demod_Rx_bit(4*c) = 0;
        elseif demod_Rx_symbol_Im(c) == 3
            demod_Rx_bit(4*c-2) = 1;
            demod_Rx_bit(4*c) = 1;
        end
    end

    biterror(ebn0index) = 0;
    for i = 1:N*k
        if demod_Rx_bit(i) ~= bitstream(i)
            biterror(ebn0index) = biterror(ebn0index) + 1;
        end
    end
    ber = biterror/(N*k);
end

theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(EbN0_dB/10))));

figure(1)
semilogy(EbN0_dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(EbN0_dB,ber,'mx-','LineWidth',2);
axis([0 15 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 16-QAM with DSSS in LOS')

figure(2)
semilogy(EbN0_dB,ser,'mx-','LineWidth',2);
axis([0 15 10^-5 1])
grid on
legend('simulation');
xlabel('Eb/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 16-QAM with DSSS in LOS')