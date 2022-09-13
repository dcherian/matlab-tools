% 
%   function [] = dchist2d(U,V,nbins)
%

function [] = dchist2d(U,V,nbins)

    merge = [U,V];
    
    vXEdge = linspace(min(U),max(U),nbins);
    vYEdge = linspace(min(V),max(V),nbins);
    mHist2d = hist2d(merge,vYEdge,vXEdge);

    nXBins = length(vXEdge);
    nYBins = length(vYEdge);
    vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
    vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
    pcolor(vXLabel, vYLabel,mHist2d); colorbar
    beautify;