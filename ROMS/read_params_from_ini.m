% [params] = read_params_from_ini(fname)
% Pass either *.in file or directory

function [params] = read_params_from_ini(fname)

    if isempty(strfind(fname,'_ini.nc'))
        fname1 = roms_find_file(fname,'ini');
        fname = [fname fname1]; % there can only be one hit
    end
    info = ncinfo(fname);
    n = length(info.Variables);
    
    for i=1:n
        name = info.Variables(i).Name;
        if strfind(name,'.')
            eval(['params.' name ' = ncread(''' fname ''',''' name ''');']);
        end
    end