% Removes main 4 tidal constituents. Returns detided signal and amplitudes
% of the tides at tide_freq.
%         [detided,error,tide_amp,tide_phase,tide_sig] = dcdetide(data,tide_freq)
%                           M2     S2     K1     O1 
% Default tide_freq = 1./[12.42, 12.00, 23.93, 25.82];
%         tide_amp  = Tidal amplitudes size(tide_freq);
%         tide_pha  = Tidal phase atan2
%         tide_sig  = Tidal signal i.e., detided + tide_sig = data
%         error     = Error in fits assuming white noise
% Assumes data timestep is 1 unit and frequency is specified in 1/unit
    
% CHANGELOG
% Added uncertainity calculations from Wunsch(1998) pg. 116     05 Dec 2011
% Modified to return fitted tide amplitude, phase               26 Nov 2011
% Now works for both column and row matrices                    23 Nov 2011

function [detided,error, tide_amp,tide_pha,tide_sig] = dcdetide(data,tide_freq)
    
    N  = length(data);
    Nt = length(tide_freq);
    
    if size(data,1) == 1
        flag1 = 1;
        data = data';
    else
        flag1 = 0;
    end
    
    if ~exist('tide_freq','var')
        tide_freq = 1./[12.42, 12.00, 23.93, 25.82];
    end
    
    A = [cos(repmat(2*pi*tide_freq,N,1) .* repmat([1:N]',1,length(tide_freq))) , sin(repmat(2*pi*tide_freq,N,1) .* repmat([1:N]',1,length(tide_freq)))];
    X = A\data;
    P = var(data)*inv(transpose(A)*A); % P = cov(x,x) : assumes white noise
    C = diag(P);
    
    % extract sin and cosine amplitudes and uncertainities
    cosX = X(1:Nt);
    sinX = X(Nt+1:Nt*2);
    cosE = C(1:Nt);
    sinE = C(Nt+1:Nt*2);

    % get results
    detided  = data - A*X;
    tide_amp = sqrt(cosX.^2 + sinX.^2);
    tide_pha = atan2(sinX,cosX);
    tide_sig = A*X;
    error    = sqrt(cosE.^2+sinE.^2);
    
    if flag1 == 1
        detided  = detided';
        tide_amp = tide_amp';
        tide_sig = tide_sig';
    end