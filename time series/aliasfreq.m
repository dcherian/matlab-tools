% returns apparent frequency given a known frequency
%       [s_app] = aliasfreq(s_true, s_sample)

function [s_app] = aliasfreq(s_true, s_sample)

    lower = s_true/s_sample - 0.5;
    upper = s_true/s_sample + 0.5;
    
    if ceil(lower) == floor(upper)
        s_app = abs(s_true - ceil(lower)*s_sample);
    else
        fprintf('\n\n problem in aliasfreq.m. Cant find an integer\n\n');
        s_app = NaN;
    end