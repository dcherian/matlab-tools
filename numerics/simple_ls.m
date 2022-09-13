% Simple Least Squares by right division
%       [x] = simple_ls(y,E)

function [x] = simple_ls(y,E)    
      x = y'/E';
%     Av = vgrid(d_range,t_range)'/Fu';
%     At = eta'/Fw';