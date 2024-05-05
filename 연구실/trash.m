clc;clear;
% 임의의 신호 생성
t = linspace(0, 1, 1000);
signal = cos(2 * pi * 50 * t) + 0.5 * cos(2 * pi * 80 * t);
% Haar 웨이블릿 변환 실행
[maxScale, step] = deal(100, 1);  % 최대 스케일과 스텝 설정
[coeffs, scales, positions] = customHaarCWT(signal, maxScale, step);
% 웨이블릿 계수 시각화
figure;
imagesc(positions, scales, abs(coeffs));
xlabel('Time');
ylabel('Scale');
title('Haar Wavelet Transform Coefficients');
colorbar;  % 색상 막대 추가