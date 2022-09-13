% finds indices along dim in array 'in' where in == val
%   function [out] = arrfind(in,val,dim)

function [out] = arrfind(in,val,dim)
    if dim == 2, in = in'; end
    
    out = bsxfun(@times,double(in == val),[1:size(in,1)]');
    
    if dim == 2, out = out'; end