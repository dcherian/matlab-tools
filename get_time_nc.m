function [time_datenum] = get_time_nc(infile)
    
    % time info
    T.ocean_time = nc_varget(infile,'ocean_time'); % time in seconds
    T.dstart = nc_varget(infile,'dstart'); % day of the first save of this run
    T.dstart_units = nc_attget(infile,'dstart','units');
    
    % NOTE T.dstart_units is a string giving the time reference
    % e.g. seconds since 2004-01-01 00:00:00
    iii = find(T.dstart_units == '-');
    iii = iii(1);
    year0 = str2num(T.dstart_units(iii-4:iii-1));
    month0 = str2num(T.dstart_units(iii+1:iii+2));
    day0 = str2num(T.dstart_units(iii+4:iii+5));
    iii = find(T.dstart_units == ':');
    iii = iii(1);
    hour0 = str2num(T.dstart_units(iii-2:iii-1));
    minute0 = str2num(T.dstart_units(iii+1:iii+2));
    second0 = str2num(T.dstart_units(iii+4:iii+5));
    datenum0 = datenum(year0,month0,day0,hour0,minute0,second0);
    time_datenum = datenum0 + T.ocean_time/86400;