% Generates series [1:n] of length of 'i'-th dimension of 'var1'

function [series] = gen_ser(var1,i)
    series = [1:size(var1,i)];