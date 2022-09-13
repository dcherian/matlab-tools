% adds string to legend
%       [legh] = addlegend(handle,str)
% from http://www.mathworks.com/support/solutions/en/data/1-181SJ/?solution=1-181SJ

function [legh] = addlegend(handle,str,loc)

    if ~exist('loc','var') || isempty(loc), loc = 'NorthEast'; end

    % get legend handles
    hleg = legend;
    if isempty(outh)
        [~,~,outh,outm] = legend('asdfg');
    end

    % remove default title
    if length(outh) == 1 && strcmp(outm{1},'asdfg')
        outh = [];
        outm{1} = str;
        legh = legend([outh;handle],outm{:},'');
        return;
    end
    
    legh = legend([outh;handle],outm{:},str);
    set(legh,'Location',loc);