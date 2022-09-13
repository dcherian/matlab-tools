% writes setup parameters to _ini.nc file for ROMS runs
%   [] = write_params_to_ini(ininame, param, varname)

function [] = write_params_to_ini(ininame, param, varname,ncid)

    if ~exist('varname','var'), varname = inputname(2); end
    
    if ~exist('ncid','var') || isempty(ncid)
        % open file for write
        ncid = netcdf.open(ininame,'WRITE');
        ncid_provided = 0;
    else
        ncid_provided = 1;
    end
    
    try
        netcdf.reDef(ncid);
    catch ME % already in define mode
    end
    % define dimensions
    try
        dimid_const = netcdf.defDim(ncid,'const_dim',1);
    catch ME
        dimid_const = netcdf.inqDimID(ncid,'const_dim');
    end
    
    if isstruct(param)
        names = fieldnames(param);
        for i=1:length(names)
            field = param.(names{i});            
            if ~isempty(field) && (numel(field) < 2 || ischar(field))
                if isstruct(field)
                    write_params_to_ini(ininame,field,[varname '.' names{i}],ncid);
                    continue;
                end
                if islogical(field), field = double(field); end
                writename = [varname '.' names{i}];
                if ischar(field)
                    try
                        dimidstr = netcdf.defDim(ncid,[writename '_dim'],length(field));
                    catch ME
                        dimidstr = netcdf.inqDimID(ncid,[writename '_dim']);
                    end
                    dimid = [dimid_const dimidstr];
                    vartype = 'char';
                else
                    if isscalar(field)
                        dimid = [dimid_const dimid_const];
                    else
                        dimid = [];
                        warning('not setup for non-scalar/non-char constants');
                    end
                    vartype = 'double';
                end
                try
                    varnum = netcdf.defVar(ncid,writename,vartype,dimid);
                catch ME % variable defined?
                    varnum = netcdf.inqVarID(ncid,writename);
                end

                netcdf.endDef(ncid);
                netcdf.putVar(ncid,varnum,field);
                netcdf.reDef(ncid);
            end
        end
    end
    
    if ~ncid_provided
        netcdf.close(ncid);
    end
