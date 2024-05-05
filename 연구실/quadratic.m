function [tau, p, quadratic] = quadratic(input, sample_shift)
    tau=-(length(input)-1):1:length(input)-1; 
    p=length(input)+sample_shift;
    if sample_shift > 0  %right shift
        shifted_input(tau(p)+1:length(input))=input(1:length(input)-tau(p));
        shifted_input(1:tau(p)) = 0;
    elseif sample_shift < 0 %left sihft
        shifted_input(1:length(input)+tau(p))=input(-tau(p)+1:length(input));
        shifted_input(length(input)+tau(p)+1:length(input)) = 0;
    else %no shift
        shifted_input(:)=input;
    end
    quadratic=input.*conj(shifted_input); %x(t)x*(t-tau)
end 
