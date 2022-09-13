% Calculate streamfunction on C-grid (interior RHO points)
%       [psi,xpsi,ypsi,zpsi] = streamfn_cgrid(grid,ubar,vbar,H)
%               - H is TOTAL water depth (zeta + bathymetry) on ALL RHO points
% supply grid structure with 
%       xu,yu,zu & xv,yv,zv & xr,yr,zr (H on xr,yr) - ALL vectors

function [psi,xpsi,ypsi] = streamfn_cgrid(grid,ubar,vbar,H)
 
    % move to RHO points & get transport
    Urho = avg1(ubar(:,2:end-1),1) .* H(2:end-1,2:end-1);
    Vrho = avg1(vbar(2:end-1,:),2) .* H(2:end-1,2:end-1);
    
    xpsi = grid.xr(2:end-1);
    ypsi = grid.yr(2:end-1);
    
    if size(xpsi,1) == 1, xpsi = xpsi'; end
    if size(ypsi,1) == 1, ypsi = ypsi'; end
    
    mask = isnan(Urho).*isnan(Vrho);
    
    % cumulative integration from right to left
    intv = fillnan(flipdim(cumtrapz(flipud(xpsi),flipdim(repnan(Urho,0),1),1),1) .* ~mask,0); 
    
    psi = intv;
    
function [um] = avgy(um)
    um = (um(:,1:end-1,:,:)+um(:,2:end,:,:))/2;

function [um] = avgx(um)
    um = (um(1:end-1,:,:,:)+um(2:end,:,:,:))/2;

function [um] = avgz(um)
    um = (um(:,:,1:end-1,:)+um(:,:,2:end,:))/2;
