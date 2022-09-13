% Finds a filename in the current directory using regexp()

function [fname] = find_file(pattern)
    
    list = ls;
    index = regexp(cellstr(list),regexptranslate('wildcard',pattern));
    fname = list(~cellfun('isempty',index),:);