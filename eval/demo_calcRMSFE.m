clear;
%-------------------------------------------------------------------------%
%- This code compares two ways to calculate the RMSFE over T periods and
%- M draws of the forecast errors. The first calculates the RMSFE for each
%- draw, the second averages the squared errors, then averages over period
%- and then takes the square root. While the former is the correct way to
%- calculate the RMSFE, the goal of this code is to show that the
%- differences between the two calculations are virtually non-existent, if
%- you look at the mean of the RMSFE across draws in the case of the first
%- approach while there are small differences when you look at the
%- percentiles. These statements are based on the assumption of Normally
%- distributed forecast errors. 
%-------------------------------------------------------------------------%

T = 100; % # of periods
M = 1000; % # of draws
A = randn(T, M) * 1 + 5; % forecast errors
A2 = A .^ 2; % squared forecast errors

% compare mean and percentiles
mean1_A = mean(sqrt(mean(A2, 1))); % mean over RMSFE across draws -> new way
mean2_A = sqrt(mean(mean(A2, 2))); % square root of the mean absolute error, averaged across draws
prctiles1_A = prctile(sqrt(mean(A2, 1)),[5, 50, 95]);
prctiles2_A = sqrt(prctile(mean(A2, 2),[5, 50, 95]));