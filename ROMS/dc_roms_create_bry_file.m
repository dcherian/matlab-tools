% Function that takes in structure, writes out appropriate grid and passive
% tracer  information to boundary conditions file.
% inspird by d_obc_roms2roms.m

function [] = dc_roms_create_bry_file(S)

	OBC.west  = S.boundary(1);
    OBC.east  = S.boundary(2);
    OBC.south = S.boundary(3);
    OBC.north = S.boundary(4);
    
    VarGrd = {'spherical',                                                ...
          'Vtransform', 'Vstretching',                                ...
          'theta_s', 'theta_b', 'Tcline', 'hc',                       ...
          's_rho', 'Cs_r', 's_w', 'Cs_w'};
    
    % write out to file
    for ii=1:length(VarGrd)
        if isfield(S,VarGrd{ii})
            nc_write(S.ncname,VarGrd{ii},S.(VarGrd{ii}));
        end
    end
    try
        nccreate(S.ncname,'bry_time','dimensions',{'bry_time' 1});
    catch ME
        warning('bry_time already exists?');
    end
    ncwriteatt(S.ncname,'bry_time','calendar','none');
    
    bry = {'west','east','south','north'};
    out_var = {'u','v','rho'};
    if S.NPT > 0
        pt_dim = {{'eta_rho' S.Mm+2 's_rho' S.N 'bry_time' 1}; ...
                  {'eta_rho' S.Mm+2 's_rho' S.N 'bry_time' 1}; ...
                  { 'xi_rho' S.Lm+2 's_rho' S.N 'bry_time' 1}; ...
                  { 'xi_rho' S.Lm+2 's_rho' S.N 'bry_time' 1}};
    end
    
    if ~S.spherical
        ax = {'x','y'};
    else
        ax = {'lon','lat'};
    end
    
    for ii=1:length(bry)
        % if needed
        if OBC.(bry{ii})
            % write out grid for OBC
            for jj=1:length(out_var)
                for kk=1:length(ax)
                    part_name = [ax{kk} '_' out_var{jj}];
                    full_name = [ax{kk} '_' out_var{jj} '_' bry{ii}];

                    switch(ii)
                        case 1 % west
                            var = S.(part_name)(1,:);
                        case 2 % east
                            var = S.(part_name)(end,:);
                        case 3 % south
                            var = S.(part_name)(:,1);
                        case 4 % north  
                            var = S.(part_name)(:,end);
                    end               

                    if size(var,1) == 1, var = var'; end

                    nc_write(S.ncname,full_name,var);     
                end
            end
            
            % write out passive tracer for OBC
            for mm=1:S.NPT
                varname = sprintf('dye_%s_%02d',bry{ii},mm);
                nccreate(S.ncname,varname,'dimensions',pt_dim{ii});
                
                % maybe add more attributes
                ncwriteatt(S.ncname,varname,'time','bry_time');
                ncwriteatt(S.ncname,varname,'units','kilogram meter-3');
            end
        end                
    end   
    