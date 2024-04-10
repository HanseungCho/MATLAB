a=rand(10,1);
x=sum(a,1);
[N,edges]=histcounts(x);

tmp1=edges(1:end-1);
tmp2=edges(2:end);
xaxis =(tmp1+tmp2)/2;

plot(xaxis,N)

xmean = mean(x);
xstd = std(x);

z = (x-xmean)/xstd; 
[N,edges]=histcounts(x);
tmp1=edges(1:end-1);s
tmp2=edges(2:end);
xaxis=(tmp1+tmp2)/2;

plot(xaxis,N)

prob = N/10^7;
probdensity = prob/(edges(2)-edges(1));

plot(xaxis,probdensity)

theo = normpdf(xaxis,0,1)
