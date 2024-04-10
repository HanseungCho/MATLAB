clc;
clear;
load('44.1sound(ncs).mat');
x=signal(1:4:end);

v=max(abs(x));
out=compand(x,255,v,'mu/compressor');
codebook=linspace(min(out),max(out),16);
diff=codebook(2)-codebook(1);

for i= 1:15
    p(i)=min(out)+(2*i-1)*(diff/2);
end

[index,quants]=quantiz(out,p,codebook);

out_c=compand(quants,255,v,'mu/expander');
sound(signal,44100);
sound(out_c,11025)
subplot(2,1,1);
plot(x);
subplot(2,1,2);
plot(out_c);       






