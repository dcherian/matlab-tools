function [result] = check_gap(data, ind, ln)

    stop_ind = ind + repmat(ln,length(ind),1);
    s = size(data);
    
    for i =1:s(2)
        for j=1:length(ln)
            for k=1:length(ind)
                if ind(k,j) ~= 0
                    nan = isnan(data(ind:stop_ind,i));
                    if find(nan == 1)
                        fprintf('\n Invalid Indices! \n');
                        result = 0;
                        break
                    end
                end
            end
        end
    end
    
    fprintf('\n Indices look alright.');
    result = 1;