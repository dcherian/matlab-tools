% [xax,yax,zax,xunits,yunits] = gcm_var_grid(fname,varname)

function [xax,yax,zax,xunits,yunits] = gcm_var_grid(fname,varname)

    % determine where on the C-grid these values lie 
    varcoords = nc_attget(fname,varname,'coordinates');
    if ~isempty(findstr(varcoords,'XU'))
      pos = 'u';
    elseif ~isempty(findstr(varcoords,'XV'))
      pos = 'v';
    elseif ~isempty(findstr(varcoords,'XC'))
      pos = 'c';
    else
      error('Unable to parse the coordinates variables to know where the data fall on C-grid')
    end
    
    if strcmp(varname,'W'), pos = 'w'; end

    switch pos
        case 'u'
            xax = ncread(fname,'Xp1');
            yax = ncread(fname,'Y');
            zax = ncread(fname,'Z');
            
            xunits = ncreadatt(fname,'Xp1','units');
            yunits = ncreadatt(fname,'Y','units'); 
            
        case 'v'
            xax = ncread(fname,'X');
            yax = ncread(fname,'Yp1');
            zax = ncread(fname,'Z');
            
            xunits = ncreadatt(fname,'X','units');
            yunits = ncreadatt(fname,'Yp1','units');
            
        case 'w'
            xax = ncread(fname,'X');
            yax = ncread(fname,'Y');
            zax = ncread(fname,'Zl');
            
            xunits = ncreadatt(fname,'X','units');
            yunits = ncreadatt(fname,'Y','units'); 

        otherwise            
            xax = ncread(fname,'X');
            yax = ncread(fname,'Y');
            
            if strcmp(varname, 'Eta');
                zax = [];
            else
                zax = ncread(fname,'Z');
            end
            
            xunits = ncreadatt(fname,'X','units');
            yunits = ncreadatt(fname,'Y','units'); 
    end