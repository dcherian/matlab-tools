function [iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,tend)
    if isempty(tindices)
        tindices = [1 Inf];
        dt = 1;
    else
        switch length(tindices)
            case 1
                dt = 1;
                tindices(2) = tindices(1);

            case 2
                dt = 1;

            case 3
                dt = tindices(2);
                tindices(2) = tindices(3);
                tindices(3) = [];
        end
    end

    if tindices(2) < tindices(1)
        tindices(2) = tindices(1)+tindices(2);
    end
    if isinf(tindices(2)), tindices(2) = tend; end
    if isinf(tindices(1)), tindices(1) = tend; end

    stride = [1 1 1 dt];

    if (tindices(2)-tindices(1)) == 0
        iend = 1;
        dt = tindices(2);
    else
        if ~isinf(slab)
            iend   = ceil((tindices(2)-tindices(1))/slab/dt);
        else
            iend = 1;
        end
    end
    % END mod_movie section
    nt = ceil((tindices(2)-tindices(1)+1)/dt);
    
    
