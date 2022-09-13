function [] = fix_subplot2x2(h)
    if length(h) ~= 4
        fprintf(' \n need to pass 4 handles');
    end
    
    % p =  [left, bottom, width, height] 
    
    % Top left
    p = get(h(1), 'pos');
    p(1) = p(1) - 0.05;
   % p(3) = p(3) + 0.1;
    set(h(1), 'pos', p);
    
    % Top Right
    p = get(h(2), 'pos');
    p(1) = p(1) - 0.08;
    %p(3) = p(3) + 0.1;
    set(h(2), 'pos', p);
    
    % Bottom left
    p = get(h(3), 'pos');
    p(1) = p(1) - 0.05;    
    p(2) = p(2) + 0.05;
    %p(3) = p(3) + 0.1;
    set(h(3), 'pos', p);
    
    % Bottom Right
    p = get(h(4), 'pos');
    p(1) = p(1) - 0.08;
    p(2) = p(2) + 0.05;
    %p(3) = p(3) + 0.1;
    set(h(4), 'pos', p);
end