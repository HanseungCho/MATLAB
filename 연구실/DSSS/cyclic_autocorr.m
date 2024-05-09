function [cf, tau, psd, p] = cyclic_autocorr(input, sample_shift, fs)
    tau=-(length(input)-1):1:length(input)-1; 
    p=length(input)+sample_shift;
    shifted_input=circshift(input, sample_shift);
    quadratic=input.*conj(shifted_input); %x(t)x*(t-tau)
    caf=fftshift(fft(quadratic));
    cf=linspace(-fs/2, fs/2-fs/length(caf), length(caf));
    psd=abs(caf);
end 

%plot x, y axis : xlabel('x=cylic frequency') ylabel(sprintf('y=cyclic autocorrelation (tau = %d)', tau(p)*(1/fs)))