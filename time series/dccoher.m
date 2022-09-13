% Calculates and returns frequencies, coherence (complex), phase (in degrees)
%       [fc,C, theta] = dccoher(x,y,dt,n,plot)
%       n -> number of frequency bands to average.
% Note: frequency units are 1/dt, plots squared coherence and phase

function [fc,C,theta] = dccoher(x,y,dt,n,plot)

    % "To remove delta functions from the origin, the mean must be removed rom
    % the series" - Bendat & Piersol pp. 180
    x = x(:)-nanmean(x(:));
    y = y(:)-nanmean(y(:));

    [f1,~,xhat] = dcfft(x,dt);
    [f1,~,yhat] = dcfft(y,dt);  
    
    if mod(n,2) == 0
        fprintf('\n Warning: n is even. Increased by 1.\n');
        n = n+1;
    end
    
    filter = ones(n,1)/n;
    
    phi_xy = conv(yhat.*conj(xhat),filter,'same');
    phi_yy = conv(yhat.*conj(yhat),filter,'same');
    phi_xx = conv(xhat.*conj(xhat),filter,'same');
    
    f2 = conv(f1,ones(n,1)/n,'same');
    
    % coherence
    C = phi_xy(n:n:end-1) ./ sqrt(phi_xx(n:n:end-1) .* phi_yy(n:n:end-1));
    theta = angle(C)*180/pi;
    fc = f2(n:n:end-1);
    % squared coherence
    C2 = abs(C).^2;
    
    sig = coher_sig(0.05,n);
    
    if exist('plot','var')
        figure;
        subplot(2,1,1)
        semilogx(fc,C2,'LineWidth',1.5);
        hold on;
        semilogx([get(gca,'XLim')],[sig sig],'k--','LineWidth', 1.5);
        title([' Coherence Squared. ', num2str(sum(C2 > coher_sig(0.05,n))/length(C2)*100), '% values above 95%']);
        beautify;
        
        subplot(2,1,2);
        semilogx(fc,theta,'LineWidth',1.5);
        logliney(0);
        title([' Phase.']);
        beautify;
    end