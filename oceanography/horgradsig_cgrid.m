% calculates horizontal gradient of variable 'var' along axis 'ax' ON SIGMA LEVEL
% using co-ordinate transformation as in WikiROMS
% does _forward_ difference 
%    [horgrad_sig] = horgradsig_cgrid(rgrid,vgrid,var,ax1)
%         rgrid - from roms_get_grid.m - IMPORTANT - expects [z,y,x]
%             - needs rgrid.zw, rgrid.s_w
%         vgrid - has vgrid.xmat,vgrid.ymat,vgrid.zmat for var - expects [x,y,z]
%               and also appropriate vgrid.s vector (S.s_w or S.s_rho)
%         var  - variable [x y z t]
%         ax1 - 1,2 for x,y

function [horgrad_sig] = horgradsig_cgrid(rgrid,vgrid,var,ax1)

error('NOT WORKING YET');

switch ax1
    case 1
        axmat = vgrid.xmat;
        ax2 = 2;
    case 2 
        axmat = vgrid.ymat;
        ax2 = 1;
end

if size(rgrid.zw,1)-1 == size(axmat,3)
    rgrid.zw = permute(rgrid.zw,[3 2 1]);
end

% Hz = dz/d?
% We compute Hz discretely as ?z/?? since this leads to the vertical sum of Hz
% being exactly the total water depth D. (from WikiROMS / Manual)
Hz = bsxfun(@rdivide,diff(rgrid.zw,1,3),diff(permute(rgrid.s_w',[3 2 1])));

% adjust Hz for non-RHO point variables
if size(Hz,1) ~= size(var,1), Hz = avg1(Hz,1); end
if size(Hz,2) ~= size(var,2), Hz = avg1(Hz,2); end
    
% (dz/dx)_?
dzdx_s = diff(vgrid.zmat,1,ax1)./diff(axmat,1,ax1);

% (df/dx)_?; x = ax1
dfdx_z = bsxfun(@rdivide,diff(var,1,ax1),diff(axmat,1,ax1));

% df/d?
dfds = nan(size(var));
dfds(:,:,2:end-1,:) = avg1(bsxfun(@rdivide,diff(var,1,3),permute(diff(vgrid.s'),[3 2 1])),3);
% pad on bottom and surface values.
dfds(:,:,1,:) = dfds(:,:,2,:);
dfds(:,:,end,:) = dfds(:,:,end-1,:); 

% chain rule power!
horgrad_sig = dfdx_z+ bsxfun(@times,avg1(1./Hz,ax1) .* dzdx_s, avg1(dfds,ax1));
