% given distribution return
function [interval] = calc95(in)

    in = sort(in);
    interval = in([floor(0.025*length(in)) ...
                   ceil(0.975*length(in))]);