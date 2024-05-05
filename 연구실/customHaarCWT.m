function [coeffs, a, b] = customHaarCWT(signal, Tb, fs)
    a = fix(Tb/1000*fs:1:Tb*fs);  
%심볼길이는 안다고 가정 Tb=1/(10^5) 대략 processing gain이 1000에서 1까지라고 잡고 search
    coeffs=zeros(length(a), length(signal));
    for s = 1:length(a)
        scale = a(s);
        % Haar 웨이블릿 정의
        wavelet = [ones(1, scale), -ones(1, scale)];
        wavelet = wavelet / sqrt(sum(wavelet.^2));
        %coeff = conv(signal, wavelet);
        %coeffs(s,:) = coeff(length(wavelet):end-(length(wavelet)-1));  
        coeff = conv(signal, wavelet, 'same');  % 'same' 옵션으로 길이 일치
        coeffs(s,:) = coeff;  % 할당 오류 방지
    end
    b = 0:1*(1/fs):length(coeffs)*(1/fs);
end
