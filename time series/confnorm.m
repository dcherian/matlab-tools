% Confidence intervals for a normally-distributed variate
%       [lower,upper] = confnorm(alpha,N)
% For 95% confidence, alpha = 0.05
%                       N = length(data)
% Bendat & Piersol 4th Ed. pp. 101
% Use for correlations

function [lower, upper] = confnorm(alpha,N)

    beta  = 2/sqrt(N-3)*abs(norminv(alpha/2));
    lower = (exp(-beta)-1)/(exp(-beta)+1);
    upper = (exp(beta)-1)/(exp(beta)+1);