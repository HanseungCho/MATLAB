function cvv = cyclic_covarianceV(input, sample_shift, t, fixed_cf)
    rv=[];
    iv=[];
    for k=1:length(sample_shift)
        shifted_input=circshift(input, sample_shift(k)); %sample_shift는 벡터
        quadratic=input.*shifted_input; %x(t)x*(t-tau)
        caf=(1/length(input))*fft(quadratic);
        rv=[rv real(caf)];
        iv=[iv imag(caf)];      
    end
    cvv=[rv iv];
end 