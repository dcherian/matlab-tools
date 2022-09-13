function [ind,num,spillover] = find_gap(var,ln)

% [ind,num,spillover] = find_gap(var,ln)
% 
% Finds and returns indices of sequential valid data locations of lengths 'ln'
% in the rows of 'var'. Returns:
%           ind -> (max(num),length(ln)); stores start indices for range. 
%               end index may be obtained from ind(:,i) + ln(i)
%           num(i) -> number of possible choices of length ln(i);
%           spillover -> stores last valid index in examined data record length

    s = size(var);
    nan = isnan(var);

    if s(2) ~= 1
        % make a compounded isnan variable and recursively call find_gap
        nan1 = isnan(var(:,1));
        for i=2:s(2)
            nan1 = nan1 | isnan(var(:,i));
        end
        % Repopulate NaN's
        nan1 = double(nan1);
        nan1(find(nan1 == 1)) = NaN;
        [ind,num,spillover] = find_gap(nan1,ln);
        return
    end % if
    
    % start(i):stop(i) gives contiguous NaN locations.
    start = find(diff(nan) == 1) + 1;
    stop = find(diff(nan) == -1);

    % Special cases   
    s1 = size(start); s2 = size(stop);
    if isnan(var(1)) == 1 % first gap starts at the beginning of var
        start(1:end+1)=[1;start];
    end
    if isnan(var(end)) == 1 % last gap extends till end of var
        stop(end+1) = s(1);
    end
    if s1(1) == 0 & s2(1) == 0
        fprintf('\n No gaps.');
		for i=1:length(ln)
            num(i) = floor(s(1)/ln(i));
            ind(1:num(i),i) = 1 + [0:floor(s(1)/ln(i))-1]'*ln(i);
            spillover(1:num(i),i) = s(1);
        end
        return    
    end
    
    % Pre-allocate and initialize variables
    gaps = [start(1)-1;start(2:end)-stop(1:end-1);s(1)-stop(end)];
    ind = zeros(sum(floor(gaps./min(ln))),length(ln));
    num = zeros(length(ln),1);
    spillover = ind;
    
    start1 = 1;
    j = 1;
    start(end+1) = s(1);
    stop(end+1) = NaN; % dummy value
    
    for i=1:length(ln)
        for j=1:length(start)
            ss = start(j)-start1+1;
            if ss >= ln(i)
                n_prev = num(i);
                num(i) = n_prev + floor(ss/ln(i));
                ind(n_prev+1:num(i),i) = start1 + [0:floor(ss/ln(i))-1]'*ln(i);
                spillover(n_prev+1:num(i),i) = start(j)-1;
            end
            start1 = stop(j) + 1;
        end % j loop
    end % i loop
    
    % Check whether it's working
    check_gap(var,ind,ln);  