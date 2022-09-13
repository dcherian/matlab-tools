% Calculates and plots structure function

function [sf] = dcstructfn(data)
    
    data = data - mean(data);
    sf = 2*var(data).*(1-xcorr(data,data));
    
    figure
    plot(sf)
    title('Structure Function')