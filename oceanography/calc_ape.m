% pass each 2d section to this script <-- THEEK KAR!
% assumes flat bottom
%
%       [ape] = ape(xpsi,zw,rho)
%           xpsi - one horizontal axis vector (PSI points for corners)
%           zw - vertical axis vector (W points for corners)
%           dy - grid spacing along other horizontal axis
%           rho - density field = length(xr) x length(zr)

function [APE] = calc_ape(xpsi,zw,rho)
    
    rho = squeeze(rho);
    rho = rho(2:end-1,:); % interior RHO points only
    nbins = 50;
    R0 = 1023;
    
    if size(xpsi,1) == 1, xpsi = xpsi'; end
    if size(zw,1) == 1, zw = zw'; end

    bins = linspace(min(rho(:)),max(rho(:)),nbins);
    
    % initial grid
    [ximat,zimat] = meshgrid(xpsi,zw); 
    
    % RHO points
    [xrmat,zrmat] = meshgrid(avg1(xpsi),avg1(zw)); 
    xr = avg1(xpsi); zr = avg1(zw);
    
    % figure out area of each cell in the grid
    dA = (diff(ximat(1:end-1,:),1,2) .* diff(zimat(:,1:end-1),1,1))';

    % calculate area of each density class
    for i = 1:length(bins)-1
        mask = rho >= bins(i) & rho < bins(i+1);
        binvol(i) = sum(dA(mask));
    end
    
    % calculate available volume for redistribution - should be able to account for sloping bottom here
    % hmmm...
    
    % redistribute under new "stretched" grid
    % assumes no bathymetry
    X = xpsi(end) - xpsi(1);
    dznew = binvol./X;
    znew = [zw(1), zw(1) + cumsum(dznew)];
   
    % LPE = lowest PE state
    LPE = -sum(9.81.* avg1(bins,2) .* (znew(1:end-1).^2 - znew(2:end).^2)/2 ... % int rho g z dz
                ./ R0./abs(znew(1))); %./ R0 ./ cellvol
    % IPE = Initial PE
    IPE = domain_integrate2(9.81*rho.*zrmat'./R0,xr,zr)./X./abs(znew(1));
    
    % calculate ape
    APE = IPE - LPE;
    
function [out] = domain_integrate2(in,xax,yax)
    out = squeeze(trapz(xax,trapz(yax,in,2),1));
            
%% old code
% histogram idea
%     N = numel(rho);
%     
%     nz = size(rho,2);
%     
%     nnew = n./N*nz;
%     
%     nnew .* xout;

% interpolate to new grid + histogram idea
%     if ~exist('nx','var'), nx = 100; end % number of horizontal points in regular grid
%     if ~exist('nz','var'), nz = 40;  end % number of vertical points in regular grid
%     if ~exist('nbins','var'), nbins = 50; end  % number of density classes
%     
%     znew = linspace(min(zr(:)),max(zr(:)),nz);
%     xnew = linspace(min(xr(:)),max(xr(:)),nx);
%     
%     % x,z initial
%     [xi,zi] = meshgrid(xr,zr);
%     % x,z new
%     [xn,zn] = meshgrid(xnew,znew);
%     
%     % interpolate to regular grid
%     rhonew = interp2(xi,zi,rho,xn,zn);
    