% Parse list of commands and returns a flag array with 1's if command is found.
% Separate commands in choices by ;
% Returns choices with recognized commands removed
%       [flag,choices] = parse_commands(cmd_list,choices)

function [flag,choices] = parse_commands(cmd_list,choices)
     
    flag = zeros([1 length(cmd_list)]);
    for i = 1:length(cmd_list)
        loc = strfind(choices,cmd_list{i});
        if ~isempty(loc)
            flag(i) = 1;
            choices = [choices(1:loc-1) choices(loc+length(cmd_list{i}):end)];
        end
    end