% [A] = avg1(X,index)

function [A] = avg1(X,index)
    
    if ~exist('index','var'), index = 1; end
%     if size(X,1) == 1
%         X = X';
%         m = 1;
%     else m = 0;
%     end
%     
%     A = (X(1:end-1,:) + X(2:end,:))/2;
%     if m, A = A'; end
%%
     dim = length(size(X));
     if size(X,index) == 1
         if isvector(X) && ~isscalar(X)
             A = avg1(X',index)';
         else
            A = X; 
         end
     else
         cnv = ones(1,dim);
         cnv(index) = 2;
         cnv = ones(cnv);

         A = convn(X,cnv./2,'valid');
     end
