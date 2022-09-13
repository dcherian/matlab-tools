% Correct overlapping tick labels
% Give it an axis and which tick index to remove
% also control number of decimal places using
% C style format string
%
%        correct_ticks(axis_name, format_string, tick_indices, hax)
function [] = correct_ticks(ax, format, ind, hax)

    if ~exist('ind', 'var'), ind = []; end
    if ~exist('hax', 'var'), hax = gca; end

    limstr = [upper(ax) 'Lim'];
    tickstr = [upper(ax) 'Tick'];
    labelstr = [tickstr 'Label'];

    for aa = 1:length(hax)
        ticks = get(hax(aa), tickstr);
        if ~isempty(format)
            for ii = 1:length(ticks)
                ticklabels{ii} = sprintf(format, ticks(ii));
                if ticks(ii) == 0
                    ticklabels{ii} = '0';
                end
            end
            set(hax(aa), labelstr, ticklabels);
        else
            ticklabels = hax.(labelstr);
        end

        if length(ticks) ~= length(ticklabels)
            ticks(ticks < hax(aa).(limstr)(1) ...
                  | ticks > hax(aa).(limstr)(2)) = [];
            hax(aa).(tickstr) = ticks;
        end

        if ~isempty(ind)
            if ~ischar(ind) && ~iscell(ind) % provided with index
                ind(ind > length(ticks)) = [];
                ticklabels{ind} = '';
            else
                if ~iscell(ind)
                    ind = cellstr(ind);
                end
                for kk=1:length(ind)
                    for tt=1:length(ticklabels)
                        if strcmpi(ticklabels{tt}, ind{kk})
                            ticklabels{tt} = '';
                        end
                    end
                end
            end
        end

        set(hax(aa), labelstr, ticklabels);
    end
end