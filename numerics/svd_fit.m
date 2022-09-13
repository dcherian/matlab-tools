% Use SVD method to solve Ex + [n] = y
%       [x] = svd_fit(y,E,factor)

function [x] = svd_fit(y,E,factor)
    
    [UU,lam,VV] = svd(E);
    
    K = find(diag(lam)>factor.*max(max(lam)));
    UK = UU(:,K);
    VK = VV(:,K);
    lamK = lam(K,K); 

    Fpinv=VK*(lamK\UK');
    x = (Fpinv*y)';
    
    Tu = UK*UK';
    imagesc(abs(Tu));
    caxis([0 1]);
    colorbar;
    pause(0.05);