% creates LTRANS initial location csv file based on ROMS float output
%   [] = ltrans_create_from_roms(fname,flt,rgrid)
% fname = output file name; flt = ROMS float file, rgrid

function [] = ltrans_create_from_roms(fname,flt,rgrid)

    fid = fopen(fname,'w');
    if fid == -1
        error('Couldn''t open file to write.');
    end
    
    roms = floats('roms',flt,rgrid);
    
    %roms.init = cut_nan(roms.init); don't think i need this anymore
    % needed if I'm running ltrans from a restarted run file
    disp(['rgrid.ocean_time(1) = ' num2str(rgrid.ocean_time(1)/86400) ' days']);
    sub = input('Enter initial time to subtract? ');
    disp(['subtracting ' num2str(sub) ' seconds']);
    roms.init(:,4) = roms.init(:,4) - sub;
    if any(roms.init(:,4)) < 0, error('subtracted time > initial deployment'); end
    fprintf(fid,'%.2f,%.2f,%.2f,%d \n',roms.init');
    fprintf('Added %d floats to %s \n',size(roms.init,1),fname);

    fclose(fid);