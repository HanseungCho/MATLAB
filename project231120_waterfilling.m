clc; clear; %find the waterfilling power allocation and the capacity
N0=1;
sum_Tx_power=[1, 5, 10, 50, 200, 1000];
channel_numbers=10;
raych=(randn(1,channel_numbers)+randn(1,channel_numbers)*1j)/sqrt(2);
raych_abs_sqr=abs(raych).^2; %|H(f)|^2 channel in exp dist
sorted_raych=sort(raych_abs_sqr);
p_bar=-N0./sorted_raych; %p_bar=p-1/lambda
lambda=zeros(1,length(sum_Tx_power));
p=zeros(length(sum_Tx_power),channel_numbers);
C=zeros(1,length(sum_Tx_power));
for i=1:length(sum_Tx_power)
    lambda(i)=1/((sum_Tx_power(i)-sum(p_bar))/channel_numbers);
    p(i,:)=(1/lambda(i))+p_bar;
    p(i,:)=max(p(i,:),0);
    C(i)=sum(log(1+(p(i,:).*sorted_raych/N0))/log(2));
    subplot(2,length(sum_Tx_power)/2,i);
    plot(sorted_raych,p(i,:),'-o');
    str1="Sum of transmission power:";
    str2=string(sum_Tx_power(i));
    str3="Capacity:";
    str4=string(C(i));
    str5=append(str1, str2);
    str6=append(str3, str4);
    legend(str5); 
    lgd = legend;
    lgd.Title.String = str6;
    xlabel('|H(f)|^2');
    ylabel('P(|H(f)|)');
    legend('Location','southeast')
end
