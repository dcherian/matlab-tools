% Removes columnwise mean from each column of var. works with nan's

function [A] = nan_detrend(var)
    A = var-repmat(nanmean(var),size(var,1),1);