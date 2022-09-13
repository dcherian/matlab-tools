% figure out slab size for reading from netcdf files.
% For energy script, it accounts for extra points being read in.
%           [slab] = roms_slab(fname, energy_flag, ntavg)

function [slab] = roms_slab(fname, energy_flag, ntavg)
    
    if ~exist('energy_flag','var'), energy_flag = 0; end
    
    varname = 'u';
    vinfo = ncinfo(fname,varname);
    s = vinfo.Size;
    
    slab = floor((216*480*40*11)/((s(1)-2)*(s(2)-2)*s(3)*s(4))) + 11 - 2;
    
    if energy_flag
        if ~exist('ntavg','var'), error('Supply ntavg for energy correction'); end
        slab = slab - ntavg + 1; 
    end