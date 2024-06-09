clc;clear;
cvm=[0 1 2 3 5 9 10 3];
ind0=find(cvm(1,1+1:length(cvm)) > 4);
ind1=ind0(ind0>5);