%  [f1,psd,ind] = dcpsd(data,dt,nf)
%          dt - time interval
%          nf - number of frequency bands to average

function [f1,psd,ind] = dcpsd(data,dt,nf,mark_peaks,factor)

    [freq,pow,cn] = dcfft(data,dt);
    
    if mod(nf,2) == 0
        fprintf('\n Warning: nf is even. Increased by 1.\n');
        nf = nf+1;
    end
    
    p1 = conv(abs(cn).^2,ones(nf,1)./nf,'same');    
    
    psd = length(data)*dt* p1(max(nf-1,1):nf:end);
    f1 = freq(max(nf-1,1):nf:end);
    
    figure;
    loglog(f1,psd,'b','LineWidth',1.5);
    title('PSD');
    
    if exist('mark_peaks','var') & mark_peaks ~= 0
        [pks,ind] = findpeaks(psd,'SORTSTR','descend');
        loglinex(freq(ind(1:mark_peaks)),factor);
        ind = ind(1:mark_peaks);
    else
        ind = NaN;
    end
    
    beautify;
    
    