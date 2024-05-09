function [cf, tau, psd, p] = cyclic_autocorr_ver1(input, sample_shift, fs)
    tau=-(length(input)-1):1:length(input)-1; 
    p=length(input)+sample_shift;
    if tau(p) > 0  %right shift
        shifted_input(tau(p)+1:length(input))=input(1:length(input)-tau(p));
        shifted_input(1:tau(p)) = 0;
    elseif tau(p) < 0 %left sihft
        shifted_input(1:length(input)+tau(p))=input(-tau(p)+1:length(input));
        shifted_input(length(input)+tau(p)+1:length(input)) = 0;
    else %no shift
        shifted_input(:)=input;
    end
    quadratic=input.*conj(shifted_input); %x(t)x*(t-tau)
    caf=fftshift(fft(quadratic));
    cf=linspace(-fs/2, fs/2-fs/length(caf), length(caf));
    psd=abs(caf).^2;
end 

%plot x, y axis : xlabel('x=cylic frequency') ylabel(sprintf('y=cyclic autocorrelation (tau = %d)', tau(p)*(1/fs)))