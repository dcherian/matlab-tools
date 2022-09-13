function [vars] = gcm_binread()

    vars.Eta = rdmds('Eta',NaN);
    vars.U = rdmds('U',NaN);
    vars.V = rdmds('V',NaN);
    vars.W = rdmds('W',NaN);
    vars.Temp = rdmds('T',NaN);
    vars.S = rdmds('S',NaN);
    vars.X = rdmds('XC');
    vars.Y = rdmds('YC');
    vars.Z = squeeze(rdmds('RC'));
    vars.phl = rdmds('PHL',NaN);
    vars.X = vars.X(:,1);
    vars.Y = vars.Y(1,:);