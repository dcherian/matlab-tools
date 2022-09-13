% function [out,xax,yax,zax,grd] = dc_roms_read_data(folder,varname,tindices, ...
%                        volume,stride,grd, ftype, dtype)

function [out,xax,yax,zax,grd] = dc_roms_read_data(folder,varname,tindices, ...
                        volume,stride,grd, ftype, dtype)

    disp(['Reading ' varname]);

    out = []; xax = []; yax = []; zax = [];
    % set flags
    k = 1;
    quitflag = 0;
    objflag = 0;

    % set inputs
    if ~exist('tindices','var') || isempty(tindices), tindices = [1 Inf]; end
    if ~exist('volume','var') || isempty(volume), volume = {}; end
    if ~exist('stride','var') || isempty(stride), stride = [1 1 1 1]; end
    if ~exist('grd','var'), grd = []; end
    if ~exist('ftype', 'var') || isempty(ftype), ftype = 'avg'; end
    if ~exist('dtype', 'var') || isempty(dtype), dtype = 'double'; end

    if length(tindices) == 1, tindices(2) = tindices(1); end

    if isobject(folder)
        run = folder;
        folder = run.dir;
        objflag = 1;
    end

    ticstart = tic;
    if strcmpi(varname,'pv') || strcmpi(varname, 'rv')
        files{1}='ocean_vor.nc';
    else
        % get all history files
        if isdir(folder)
            files = roms_find_file(folder, ftype);
        else
            files = {folder};
        end
    end

    %unused params
    slab = Inf;

    for ii=1:length(files)
        if isdir(folder)
            fname = [folder '/' char(files(ii))];
        else
            fname = char(folder);
        end
        if strcmpi(varname,'pv') || strcmpi(varname, 'rv')
            grd = fname;
        end

        % variable information
        vinfo  = ncinfo(fname,varname);
        nt     = vinfo.Size(end);
        ndim = length(vinfo.Size);

        if ndim == 2 & ~strcmpi(ftype, 'flt')
            nt = 1;
        end

        % float variables
        if strcmpi(ftype, 'flt')
            vol = [];
            xax = [];
            yax = [];
            zax = [];

            if ndim == 2
                stride = [1 1];
            else
                stride = 1;
            end
        else
            % all other variables
            if ii == 1
                if ndim == 3
                    stride(4) = [];
                end
                if ndim == 2
                    stride(3:4) = [];
                end
                if ndim ~= 1
                    if isempty(grd)
                        if objflag
                            grd = run.rgrid;
                        else
                            grd = roms_get_grid(fname,fname,0,1);
                        end
                    end
                    [xax,yax,zax,vol] = dc_roms_extract(grd,varname,volume,1);
                    xax = squeeze(xax); yax = squeeze(yax); zax = squeeze(zax);
                else
                    % for 1D time series data, none of this is applicable
                    vol = [];
                    xax = [];
                    yax = [];
                    zax = [];

                    stride(2:end) = [];
                end
                %[~,~,~,time,xunits,yunits] = dc_roms_var_grid(grd,varname);
            end
        end

        % process tindices (input) according to file
        [~,tnew,dt,~,tstride] = roms_tindices(tindices,slab,nt);

        % if tindices(2) is Inf and there are multiple files the case 2
        % does not execute since tnew is set by roms_tindices to be nt -
        % which comes from 1 file only.
        if isinf(tindices(end)), tnew(2) = Inf; end
        stride(end) = tstride(end);

        % Case 1 : if requested data not in this file, skip to next
        if tnew(1) > vinfo.Size(end)
            tindices = tindices - nt;
            % requested data was not in supplied file.
            if ~isdir(folder)
                error(['Data not in file: ' folder]);
            end
            continue;
        else
            % Case 2 : requested data spans 2 files
            if tnew(1) <= nt && tnew(2) > nt
                % change ending for current read
                tnew(2) = nt;
                % set next read to start from beginning of new file
                tindices(1) = 1;
                tindices(end) = tindices(end)-nt;
            else
                % Case 3 : requested data finishes in current file
                if tnew(2) <= nt && ~isinf(tindices(end))
                    quitflag = 1;
                end
            end
        end
        [start,count] = roms_ncread_params(ndim,0,1,slab,tnew,dt,vol);

        if strcmpi(dtype, 'single')
            temp = single(ncread(fname,varname,start,count,stride));
        else
            temp = double(ncread(fname,varname,start,count,stride));
        end
        if count(end) == 1 && ii == 1 && quitflag == 1
            % first file has the single timestep
            out = squeeze(temp);
            return;
        end

        % make sure appending timestep always works
        dimsave = length(start);

        if isvector(temp)
            out(k:k+length(temp)-1) = temp;
            if k == 1 & ~isequal(size(out), size(temp))
                out = out';
            end
            k = k+length(temp);
        else
            % call ndims instead of using ndim because i could be
            % extracting a slice
            switch dimsave
                case 2
                    out(:,k:k+size(temp,2)-1) = temp;
                case 3
                    out(:,:,k:k+size(temp,3)-1) = temp;
                case 4
                    out(:,:,:,k:k+size(temp,4)-1) = temp;
            end
            if count(end) ~= 1
                k = k + size(temp,dimsave);
            else
                k = k + 1;
            end
            clear temp;
        end
        if quitflag, break; end
    end
    out = squeeze(out);
    toc(ticstart);