% Wrapper script that calls roms_movie for only 1 time instant.
%       [] = roms_section(fname, varname, tindex, axis, index)

function [] = roms_section(fname, varname, tindex, axis, index)

    roms_movie(fname,varname,[tindex 1],axis,index);