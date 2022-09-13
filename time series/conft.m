% Confidence intervals for a t-distributed variate
%       [lower,upper] = conft(alpha,nf,test)
% For 95% confidence, alpha = 0.05
% Bendat & Piersol 4th Ed. pp. 90

% Test

% mean(data) + [lower,upper]*std(data)/sqrt(N);

function [lower,upper] = conft(alpha,nf,test)

    if ~exist('test', 'var'), test = 0; end

    if test
        test_conft;
        return;
    end

    lower =  tinv(alpha/2,nf);
    upper =  - tinv(alpha/2,nf);
end

function [] = test_conft

    A = [60 61 47 56 61 63 65 69 54 59 43 61 55 61 56 48 67 65 60 58 57 62 57 58 53 59 58 61 67 62 54];
    alpha = 0.1;
    nf = length(A);

    [lower,upper] = conft(alpha, nf);

    lower = mean(A) + lower * std(A)/sqrt(nf);
    upper = mean(A) + upper * std(A)/sqrt(nf);

    assert((lower - 56.85) < 1e-2)
    assert((upper - 60.37) < 1e-2)

    disp('Test passed');
end