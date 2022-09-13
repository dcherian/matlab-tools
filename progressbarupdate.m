% progressbarupdate(cpb,value,text)
% VALUE in percentage.
% TEXT can be empty or optional

function [] = progressbarupdate(cpb,value,text)

    cpb.setValue(value);
    if ~exist('text','var') || isempty(text)
        return
    else
        cpb.setText(text);
    end