function [out] = breakapart(in)

    branchPoints = bwmorph(in,'branchpoints');
    branchPoints = imdilate(branchPoints,strel('disk',3));
    % now i've broken skeleton into branches
    out = in & ~branchPoints;