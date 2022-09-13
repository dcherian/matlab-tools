% Confidence intervals for a chi-squared-distributed variate
%       [lower,upper] = confchi2(alpha,nf)
% For 95% confidence, alpha = 0.05
% Bendat & Piersol 4th Ed. pp. 90

% Test
% A = [60 61 47 56 61 63 65 69 54 59 43 61 55 61 56 48 67 65 60 58 57 62 57 58 53 59 58 61 67 62 54];
% alpha = 0.1, lower = 22.91, upper = 54.22
%                        .92            .2569
% var(data)*[lower,upper]

function [lower,upper] = confchi2(alpha,nf)
    
    upper = nf/chi2inv(alpha/2,nf);
    lower = nf/chi2inv(1-alpha/2,nf);