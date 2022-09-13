function [data] = dc_roms_zslice_var(data,depth,grd)
% $Id: roms_zslice_var.m 358 2008-04-07 14:15:03Z zhang $
% Get a constant-z slice out of a 4-D ROMS variable 
% [data,x,y] = roms_zslice_var(data,depth,grd)
% grd has grd.xax, grd.yax, grd.zax
% these and data are (x,y,z,t)
%
% John Wilkin
% modified - Deepak Cherian

data = permute(data,[3 2 1]);
z    = permute(grd.zax,[3 2 1]);

depth = -abs(depth);

% interpolate to requested depth 

% make 2d 'ribbon' out of the data
[N,L,M] = size(data);
data = data(:,1:L*M);

z = reshape(z,[size(z,1) size(z,2)*size(z,3)]);

% code lifted from omviz/scrum_zslice:

% pad out bottom and surface z values with -Inf and 0 respectively
z = [-Inf*ones([1 L*M]); z; zeros([1 L*M])];

% pad out bottom and surface data values
data = [NaN*ones([1 L*M]); data; data(N,:)];

z = flipud(z);
data = flipud(data);

% Find the indices of data values that have just greater depth than
% depth

zg_ind = find(diff(z<depth)~=0);
zg_ind = zg_ind + [0:1:length(zg_ind)-1]';
data_greater_z = data(zg_ind);
depth_greater_z = z(zg_ind);
        
% Find the indices of the data values that have just lesser depth
% than depth
zl_ind = find(diff(z>depth)~=0);
zl_ind = zl_ind + [1:1:length(zg_ind)]';
data_lesser_z = data(zl_ind);
depth_lesser_z = z(zl_ind);
        
% Interpolate between the data values.
alpha = (depth-depth_greater_z)./(depth_lesser_z-depth_greater_z);
data_at_depth = (data_lesser_z.*alpha)+(data_greater_z.*(1-alpha));
data = reshape(data_at_depth,[L M]);

data = data';

% Apply mask to catch shallow water values where the z interpolation does
% not create NaNs in the data
%mask(mask == 0) = NaN;
%data = data.*mask;
