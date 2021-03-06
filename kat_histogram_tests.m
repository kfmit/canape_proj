close all
clear all
clc

% first create two random data sets of size 1000
x1=randn(1000);     % make data 1
x2=randn(1000)+10; % make data 2 and make it "SOME OFFSET" higher than the other one
db=0.1;             % establish distance between bins
bins=-10:db:10;     % create bin edges for the histogram/ we have 200 bins, 201 edgesmac

% Now create histcounts of the two variables
% histcounts functions: val1 = the number of instances in this bin
% bin1= the bin edges, not that this comes back as the same input as bins!

% i chose to get ride of the probability normalization
[val1, bin1]=histcounts(x1,bins,'Normalization','probability');   %,,'Normalization','pdf');  
bin1_c=bin1(1:end-1)+db/2;  % cut off the extra edge so we can use it for plotting

[val2, bin2]=histcounts(x2,bins,'Normalization','probability');    %,'Normalization','pdf');  
bin2_c=bin2(1:end-1)+db/2;  % cut off the extra edge so we can use it for plotting

% make a figure that plots the raw data straight up
figure
plot(bin1_c, val1) % bin edges, number in each
hold on
plot(bin2_c, val2) % bin edges, number in each
title('The actual raw data')
xlabel('Edges of the bins')
ylabel('count in each bin')

% find the maxima of each count - MODE %%
[max1, ind_bin1]=max(val1); % max is the numerical max, ind_bin is the index of that maxima
[max2, ind_bin2]=max(val2); % ^same


% maximum distance between modes: distance between the modes
% our pairwise dist?
% this should be pretty damn close to what we set as the offset above
dist_mode=bin2_c(ind_bin2)-bin1_c(ind_bin1)

% total variation distance - literally how much distance is there between
% every count of every bin, larger number is more differnet
dist_tv=sum(abs(val1-val2))/2

% This is a figure that plots the bin edge of the maximum value and the max
% value
figure
stem(bin2(ind_bin1),max1)
hold on
stem(bin2(ind_bin2),max2)
hold off
title('plots showing distance between datasets mode')
ylabel('bin edge')
xlabel('max count of that distribution')

%%%%%%%%
% figure that plots random data and gets a real mode bc theyre integers 
figure
y = randn(1000,1);
edges = -6:.25:6;
bin = discretize(y,edges);
m = mode(bin);
edges([m, m+1]);
histogram(y,edges);
title('hist of straight up randomly generated data')