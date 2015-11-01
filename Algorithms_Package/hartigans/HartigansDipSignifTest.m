%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Calculates Hartigan's dip statistic and its significance for the 
% empirical pdf  XPDF (vector of sample values)
% This routine calls the routine 'HartigansDipTest' that actually calculates 
% the dip value. nboot is the user-supplied sample size of boot-strap.
%------------
% Written by F. Mechler (27 August 2002)
% Minor changes by Argyris Kalogeratos, 2012-2013.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dip, p_value, xlow,xup] = HartigansDipSignifTest (xpdf, nboot)

% calculate the DIP statistic from the empirical pdf
[dip,xlow,xup, ~, ~, ~, ~, ~] = HartigansDipTest(xpdf);
N = length(xpdf);

% calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
boot_dip = zeros(nboot, 1);
for i=1:nboot
   unifpdfboot = sort(unifrnd(0,1,1,N));
   unif_dip    = HartigansDipTest(unifpdfboot);
   boot_dip(i) = unif_dip;
end

boot_dip = sort(boot_dip);
p_value  = sum(dip<boot_dip) / nboot;

end


