% Load chindi variables 
%       [misc] = roms_load_misc(fname)

function [misc] = roms_load_misc(fname)

    if isdir(fname)
        dirname = fname;
        fname = roms_find_file(dirname, 'his');
        fname = [dirname '/' fname{1}];
    end

    if ~exist('fname','var'), fname = 'ocean_his.nc'; end
    
    vars = {'nl_tnu4', 'nl_visc4', 'nl_tnu2','nl_visc2','rdrg','rdrg2', ...
            'rho0','R0','Tcoef','Scoef','h','f'};
    
    for i=1:length(vars)
        try
            misc.(vars{i}) = ncread(fname,vars{i}); %eval(command);
        catch ME
            misc.(vars{i}) = NaN;
        end
    end
    