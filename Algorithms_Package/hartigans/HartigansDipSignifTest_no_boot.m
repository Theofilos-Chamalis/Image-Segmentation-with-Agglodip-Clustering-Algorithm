%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Calculates Hartigan's dip statistic and its significance for the 
% empirical pdf  XPDF (vector of sample values)
% This routine calls the routine 'HartigansDipTest' that actually calculates 
% the dip value. nboot is the user-supplied sample size of boot-strap.
% boot_dip is a bootstrap sample of size NBOOT of the dip statistic for a 
% uniform pdf of sample size N (the same as empirical pdf) which has been
% calculated in test_projected_unimodality.
%------------
% Written by F. Mechler (27 August 2002)
% Minor changes by Chamalis Theofilos, 2014-2015.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dip, p_value, xlow,xup] = HartigansDipSignifTest_no_boot (xpdf, nboot,boot_dip)

    % calculate the DIP statistic from the empirical pdf
    [dip,xlow,xup, ~, ~, ~, ~, ~] = HartigansDipTest(xpdf);
    p_value  = sum(dip<boot_dip) / nboot;
    
end


