% calculates horizontal gradient of variable 'var' along axis 'ax'
% does _forward_ difference but is setup for central difference
%   grid - has grid.xmat,grid.ymat,grid.zmat
%   var  - variable
%   ax1 - 1,2 for x,y

% MUST ACCOUNT FOR VARYING VARIABLE BUT INVARIANT GRID IN ONE DIRECTION
% NAN PROBLEM

function [horgrad] = horgrad_cgrid_interp(grid,var,ax1)

    % ax1 = axis along which derivative is being taken
    % ax2 = other axis

    % #1 - loop through every point, linearly interpolate neighbouring
    % profiles to same depth and calculate derivative
    s = size(var);
    % size for horgrad variable
    snew = s;
    snew(ax1) = s(ax1)-2;
    
    if ax1 == 1
        axmat = grid.xmat;
        ax2 = 2;
    else
        axmat = grid.ymat;
        ax2 = 1;  
    end

    % if invariant in one direction, shortcut
    flags.gridinv1 = check_invariant(grid.zmat,ax1);
    flags.gridinv2 = check_invariant(grid.zmat,ax2);
    flags.varinv1  = check_invariant( var,ax1);
    flags.varinv2  = check_invariant( var,ax2);
    
    % loop indices & set horgrad size
    ilow = 2; ihigh = s(1)-1;
    jlow = 2; jhigh = s(2)-1;    
    
    if flags.varinv2 && flags.gridinv2
        if ax1 == 1
            jlow = 2; jhigh = 2;
            repvec = [1 s(2) 1];
            horgrad = nan([s(1)-2 1 s(3)]);
        else
            ilow = 2; ihigh = 2;
            repvec = [s(1) 1 1];
            horgrad = nan([1 s(2)-2 s(3)]);
        end
    else
        horgrad = nan(snew);
    end

    % actually calculate derivative
    tic;
    for ii = ilow:ihigh
        for jj = jlow:jhigh
            % appropriate indices
            if ax1 == 1
                im1 = ii-1; ip1 = ii+1;
                jm1 = jj;   jp1 = jj;
            else
                im1 = ii;   ip1 = ii;
                jm1 = jj-1; jp1 = jj+1;
            end

            % figure out appropriate z grids
            zz   = grid.zmat(ii,jj,:);
            %zzm1 = grid.zmat(im1,jm1,:);
            zzp1 = grid.zmat(ip1,jp1,:);
            
            for kk = 1:s(3)
                % var @ minus/plus 1 point interpolated
                %vm1i = interp1(zzm1,squeeze(var(im1,jm1,:)),zz);
                vp1i = interp1(zzp1,squeeze(var(ip1,jp1,:)),zz);
                dx   = (axmat(ip1,jp1,kk) - axmat(im1,jm1,kk))/2;
                %horgrad(ii-1,jj-1,kk) = (vp1i(kk)-vm1i(kk))./dx;
                %if isnan(vm1i(kk)) % minus 1 point does not exist?
                % forward difference
                horgrad(ii-1,jj-1,kk) = (vp1i(kk)-var(ii,jj,kk))./dx;
                %end
                % NEED TO CORRECT FOR POINTS ABOVE AND BELOW GRID
            end
        end % jj loop
    end % ii loop
    toc;

    % appopriate repmat based on axis
    if flags.varinv2 && flags.gridinv2
        horgrad = repmat(horgrad,repvec);
    end

% check if grid is invariant in ax2 direction
function [flag] = check_invariant(var,ax2)
    dv  = diff(var,1,ax2);
    if mean(dv(:)) <= 1e-3
        % approximately invariant in the other dirn.
        flag = 1;
    end


% #2 - interpolate to uniform grid, differentiate and then interpolate back
% wont work because i'll need a lot of points to preserve resolution in
% shallow water
% zi = linspace(max(grid.zmat(:)),min(grid.zmat(:)),100);
% [ximat,yimat,zimat] = ndgrid(grid.xmat(:,1,1),grid.ymat(1,:,1),zi);
% 
% vari = interpn(grid.xmat,grid.ymat,grid.zmat,var,ximat,yimat,zimat);