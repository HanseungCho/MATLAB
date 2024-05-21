function cvv = F(input, sample_shift, t, fixed_cf)
    shifted_input=circshift(input, sample_shift); %sample_shift는 벡터
    quadratic=input.*shifted_input; %x(t)x*(t-tau)
    cvv=sum(quadratic.*exp(-1j*2*pi.*fixed_cf.*t));         
end 