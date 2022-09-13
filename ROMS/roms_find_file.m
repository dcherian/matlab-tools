% finds files for type = ini / bry / grd / fpos / flt / his / avg
%   [fname] = roms_find_file(dir,type)
% can return a list for last two types. fname does NOT contain input
% directory

function [fname] = roms_find_file(dirin,type)

    if ~isdir(dirin)
        if strcmpi(type,'his') || strcmpi(type,'avg')
            fname = dirin;
            return;
        else
            index = strfind(dirin,'/');
            dirin = dirin(1:index(end));
        end
    end

    if isempty(strfind(dirin,'config'))
        fname = [dirin '/config/'];
    else
        fname = dirin;
    end
    % ls gives different results on windows and linux
    %in = ls([fname '/*.in']);
    in = dir([fname '/*.in']);

    if size(in,1) > 1
        ii = 1;
        while size(in,1) > 1
            if strcmpi(in(ii,:).name,'floats.in') | ...
                    strcmpi(in(ii,:).name,'stations.in') | ...
                    ( isempty(strfind(in(ii,:).name,'rst') | ...
                              strfind(in(ii,:).name, '.#')) &  ...
                      strfind(in(ii,:).name,'.in'))
                in(ii,:) = [];
                continue
            end
            ii = ii+1;
        end
    end

    try
        in = in.name;
    catch ME
        error('Cannot find *.in file. check config folder');
    end

    % files from *.in
    if strcmpi(type,'ini') | strcmpi(type,'bry') |  strcmpi(type,'grd') ...
            | strcmpi(type, 'fpos')
        fname = ['/config/' grep_in([fname in],type)];
        return;
    end

    % floats
    if strcmpi(type,'flt')
        fnames = dir([dirin '/*_flt*.nc*']);
    end

    if strcmpi(type,'his')
        fnames = dir([dirin '/*_his*.nc*']);
        if isempty(fnames)
            fnames = dir([dirin '/*_avg*.nc*']);
            disp('Using avg files instead.');
        end
    end

    if strcmpi(type,'avg')
        fnames = dir([dirin '/*_avg*.nc*']);
        if isempty(fnames)
            fnames = dir([dirin '/*_his*.nc*']);
            disp('Using his files instead.');
        end
    end

    % convert from struct to names
    clear fname
    for kk=1:size(fnames)
        fname{kk}= fnames(kk,:).name;
    end


% runs grep on input file
function [str] = grep_in(fname,type)
    [~,p] = grep('-s', [upper(type) 'NAME == '],fname);
    if isempty(p.match) % catch FPOSNAM
        [~,p] = grep('-s', [upper(type) 'NAM = '],fname);
        % line in p.match must be processed to extract *.nc name
        str = sscanf(char(p.match),sprintf(' %sNAM = %%s', ...
                                           upper(type)));
    else
        % line in p.match must be processed to extract *.nc name
        str = sscanf(char(p.match),sprintf(' %sNAME == %%s', ...
                                           upper(type)));
    end
