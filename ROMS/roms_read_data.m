function [out] = roms_read_data(folder,varname,start,count,stride)
    
    disp(['Reading ' varname]);
    tic;
    if strcmpi(varname,'vor')
        files{1}='ocean_vor.nc';
    else
        % get all history files
        if isdir(folder)
            files = roms_find_file(folder,'his');
        else
            files = folder;
        end
    end
    k = 1;
    
    for ii=1:length(files)
        if isdir(folder)
            fname = [folder '/' char(files(ii))];
        else
            fname = folder;
        end
        if ii == 1
            vinfo  = ncinfo(fname,varname);
            dim = length(vinfo.Size);
            

            if ~exist('start','var') || isempty(start), start = ones([1 dim]); end
            if ~exist('count','var') || isempty(count), count = inf([1 dim]); end
            if ~exist('stride','var') || isempty(stride), stride = ones([1 dim]); end
        end
        temp = squeeze(double(ncread(fname,varname,start,count,stride)));
        if count(end) == 1
            out = temp;
            return;
        end
        
        if isvector(temp)
            out(k:k+length(temp)-1) = temp;
            k = k+length(temp);
        else        
            switch ndims(temp)
                case 2
                    out(:,k:k+length(temp)-1) = temp;
                    k = k+length(temp);
                case 3
                    out(:,:,k:k+size(temp,3)-1) = temp;
                    k = k+size(temp,3);
                case 4
                    out(:,:,:,k:k+size(temp,4)-1) = temp;
                    k = k+size(temp,4);
            end  
        end
    end 
    toc;
