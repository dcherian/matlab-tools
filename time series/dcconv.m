% Performs discrete convolution
%      [h] = dcconv(f,g)
% NOTE: If using to filter, make sure window is normalized!
% f and g must be of same dimensions (and 1D)

function [h] = dcconv(f,g)

    if length(f) < length(g)
        ff = f;
        f = g;
        g = ff;
    end
       
    nf = length(f);
    ng = length(g);        
    
    % reflect g
    g = fliplr(g);
    
    for i=1:nf
        gg = g(max(1,ng-i+1):ng); 
        if i < ng+1
            h(i) = sum(f(1:length(gg)) .* gg);
        else
            h(i) = sum(f(i-ng+1:i-ng+length(gg)) .* gg);
        end
    end
    
    for i=nf+1:nf+ng-1
        gg = g(1:min(ng,ng+nf-i));
        h(i) = sum(f(i-ng+1:nf) .* gg);
    end