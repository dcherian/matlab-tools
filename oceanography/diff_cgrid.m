% calculates derivatives : either vertical or
% horizontal gradient of variable 'var' along axis 'ax' ON Z LEVEL
% using co-ordinate transformation as in WikiROMS
% does _forward_ difference 
%           [der] = diff_cgrid(grid,var,ax1)
%
%   grid - has grid.xmat,grid.ymat,grid.zmat for var - expects [x,y,z]
%           and also appropriate grid.s vector (S.s_w or S.s_rho)
%           grid.zw, grid.s_w
%   var  - variable [x y z t]
%   ax1 - 1,2,3 for x,y,z

function [der] = diff_cgrid(grid,var,ax1)

% dfdz = 1./Hz .* dfds
%dfdz = dz_cgrid(grid,var);
dfdz = bsxfun(@rdivide,diff(var,1,3),diff(grid.zmat,1,3));

if ax1 == 3 % d/dz only
    der = dfdz;
    return
end

dfdz = avg1(dfdz,ax1);

% for horizontal derivatives
switch ax1
    case 1
        axmat = grid.xmat;
    case 2 
        axmat = grid.ymat;
end

% (dz/dx)_s
dzdx_s = diff(grid.zmat,1,ax1)./diff(axmat,1,ax1);

% (df/dx)_s; x = ax1
dfdx_s = bsxfun(@rdivide,diff(var,1,ax1),diff(axmat,1,ax1));

% chain rule power!
if size(dfdx_s,3) ~= size(dfdz,3) % for w
   der = avg1(dfdx_s,3) - bsxfun(@times,avg1(dzdx_s,3),dfdz);
else
    der = dfdx_s - bsxfun(@times,dzdx_s, dfdz);
end

debug = 0;

if debug
    figure
    h(1) = subplot(141);
    contourf(avg1(squeeze(grid.xmat(1,:,:)),ax1)/1000,avg1(squeeze(grid.zmat(1,:,:)),ax1), ...
        squeeze(dfdx_s(1,:,:)));
    title('(df/dx)_s'); colorbar; cx = caxis;
    h(2) = subplot(142);
    contourf(avg1(avg1(squeeze(grid.xmat(1,:,:)),ax1),2)/1000, ...
        avg1(avg1(squeeze(grid.zmat(1,:,:)),ax1),2), ...
        squeeze(dfdz(1,:,:)));
    title('df/dz'); colorbar; cx = caxis;
    h(3) = subplot(143);
    contourf(avg1(squeeze(grid.xmat(1,:,:)),ax1)/1000,avg1(squeeze(grid.zmat(1,:,:)),ax1), ...
        squeeze(dzdx_s(1,:,:)));
    title('(dz/dx)_s'); colorbar; cx = caxis;
    h(4) = subplot(144);
    contourf(avg1(squeeze(grid.xmat(1,:,:)),ax1)/1000,avg1(squeeze(grid.zmat(1,:,:)),ax1), ...
        squeeze(der(1,:,:)));
    title('(df/dx)_z = (df/dx)_s - (dz/dx)_s * df/dz'); colorbar; cx = caxis;
    linkaxes(h,'xy');
end

% moved to dz_cgrid.m

% if size(grid.zw,1)-1 == size(axmat,3)
%     grid.zw = permute(grid.zw,[3 2 1]);
% end
% 
% % Hz = dz/d?
% % We compute Hz discretely as ?z/?? since this leads to the vertical sum of Hz
% % being exactly the total water depth D. (from WikiROMS / Manual)
% Hz = bsxfun(@rdivide,diff(grid.zw,1,3),diff(permute(grid.s_w',[3 2 1])));

% adjust Hz for non-RHO point variables
% if size(Hz,1) ~= size(var,1), Hz = avg1(Hz,1); end
% if size(Hz,2) ~= size(var,2), Hz = avg1(Hz,2); end
% df/d?
% dfds = nan(size(var));
% dfds(:,:,2:end-1,:) = avg1(bsxfun(@rdivide,diff(var,1,3),permute(diff(grid.s'),[3 2 1])),3);
% % pad on bottom and surface values.
% dfds(:,:,1,:)   = dfds(:,:,2,:);
% dfds(:,:,end,:) = dfds(:,:,end-1,:); 
