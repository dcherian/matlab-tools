% Read Mellor Yamada state output

function [vars] = gcm_myread(fname)

    fname = find_file(fname);
    
    vars.X = ncread(fname,'X');
    vars.Y = ncread(fname,'Y');
    vars.Z = ncread(fname,'Z');
    vars.MYviscAr = ncread(fname,'MYviscAr');
    vars.MYdiffKr = ncread(fname,'MYdiffKr');
    vars.MYhbl = ncread(fname,'MYhbl');