input=[1 2 3 4 5];
tau=-(length(input)-1):1:length(input)-1; 

quadratic = zeros(1, length(tau)); % 결과를 저장할 배열을 초기화합니다.
fs=10^9;

if tau(p) > 0  % 오른쪽으로 이동
    shifted_input = [input(tau(p)+1:end) zeros(1, tau(p))];
elseif tau(p) < 0  % 왼쪽으로 이동
    shifted_input = [zeros(1, -tau(p)) input(1:end+tau(p))];
else  % 이동이 없는 경우
    shifted_input = input;
end
quadratic = input .* conj(shifted_input); % x(t)x*(t-tau)
caf=fftshift(fft(quadratic));
cf=linspace(-fs/2, fs/2, length(caf));
psd=abs(caf).^2;
plot(cf, psd);
