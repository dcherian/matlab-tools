%  [grid1] = dc_roms_get_grid(fname)
% returns all grid in matrix form for u,v,w,rho as xumat,yumat,zumat etc.

function [grid1] = dc_roms_get_grid(fname,varname)

    [grid1.xumat,grid1.yumat,grid1.zumat,grid1.tax,grid1.xunits,grid1.yunits] ...
        = dc_roms_var_grid(fname,'u');

    vars = {'v','w','r'};
    varnames = {'v','w','rho'};

    for i=1:length(vars)
        command = sprintf(['[grid1.x%smat,grid1.y%smat,grid1.z%smat,~,~,~]' ...
            '= dc_roms_var_grid(fname,''%s'');'], ...
            vars{i},vars{i},vars{i},varnames{i});
        eval(command);
    end