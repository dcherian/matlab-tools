function out = mat2vec(in)
    out = reshape(in,prod(size(in)),1);