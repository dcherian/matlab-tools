function [Zout] = mkgrid(Zin, start)
   
    if size(Zin,1) == 1, Zin = Zin'; end
    
    Zin = [start; Zin];
    Zout(1) = start;
    for i=1:length(Zin)-1
        Zout(i+1) = 2*Zin(i+1)-Zout(i);
    end
    
    if size(Zin(2:end)) ~= size(Zout)
        Zout = Zout';
    end
    
    if avg1(Zout) == Zin(2:end)
        fprintf('\n New grid looks alright. \n');
    else
        fprintf('\n New grid does not line up! \n');
    end