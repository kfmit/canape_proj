close all
clear all
clc
x1=randn(1000);
x2=randn(1500)+0.1;
db=0.1;
bins=-10:db:10;

[val1, bin1]=histcounts(x1,bins,'Normalization','probability');  
bin1_c=bin1(1:end-1)+db/2;
[val2, bin2]=histcounts(x2,bins,'Normalization','probability');  
bin2_c=bin2(1:end-1)+db/2;

figure
plot(bin1_c, val1)
hold on
plot(bin2_c, val2)
[max1, ind_bin1]=max(val1);
[max2, ind_bin2]=max(val2);
% maximum distance between modes
dist_mode=bin2_c(ind_bin2)-bin1_c(ind_bin1)

% total variation distance
dist_tv=sum(abs(val2-val1))/2


figure
stem(bin2(ind_bin2),max1)
hold on
stem(bin2(ind_bin2),max2)
%%%%%%%%
y = randn(1000,1);
edges = -6:.25:6;
bin = discretize(y,edges);
m = mode(bin);
edges([m, m+1])
histogram(y,edges)