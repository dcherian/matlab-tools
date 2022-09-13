% Fills data in 2D matrix along specified dimension using interp1
%       [output] = fill_data(input,dim)
% a) Assumes dim = 2 by default i.e., interpolates gaps in columns along
% each row.

function [out] = fill_data(in,dim)

    % need to move along opposite dimension to that which must be interpolated
    if isempty(dim)
        dim = 2;
        s = 1;
    else
        if dim == 1
            s = 2;
        else
            s = 1;
        end
    end

    x = 1:size(in,dim);
    out = in;

    for i = 1:size(in,s);
        if s == 1
            mask = ~isnan(in(i,:));    
            out(i,~mask) = interp1(x(mask), in(i,mask), x(~mask));   
        else    
            mask = ~isnan(in(:,i));    
            out(~mask,i) = interp1(x(mask), in(mask,i), x(~mask));   
        end
    end
