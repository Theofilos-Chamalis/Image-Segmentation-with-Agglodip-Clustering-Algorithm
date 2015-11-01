%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Demonstration script for the dip statistic.
% Creates some obviously unimodal and bimodal Gaussian distributions just 
% to show what dip statistic does.
%------------
% Written by Nic Price, 2006
% Minor changes by Argyris Kalogeratos, 2012-2013.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cent1 = ones(1,9);
Cent2 = 1:1:9;
sig   = 0.5;
nboot = 500; % the number of bootstrap uniform distributions to use

for a = 1:length(Cent1)
    
    xpdf(a,:) = sort([Cent1(a)+randn(1,200) Cent2(a)+randn(1,200)]); % allocate 200 points in each
    [dip(a), p(a)] = HartigansDipSignifTest(xpdf(a,:), nboot);
    
    subplot(3,3,a)
    hist(xpdf(a,:),-2:0.25:12)    
    title(['dip=',num2str(dip(a),3), ', p=',num2str(p(a),3)])
    xlim([-2 12])
    
end
