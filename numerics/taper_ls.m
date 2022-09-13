% Smoothed Tapered LS: 
% page 56 of Wunsch - Discrete Inverse and State Estimation with W=I
% gamma2 trades smoothness against noise
%
%       [x] = taper_ls(y,E,gamma2) solves Ex + [n] = y;

function [x] = taper_ls(y,E,gamma)
    
    % Problems.
    n_modes = size(E,2);
    
%     F13(1:n_modes,1) = 0;
%     F13(1:n_modes-1,2:n_modes) = (-1)*eye(n_modes-1);
%     F11 = eye(n_modes) + F13;
%     
%     S = gamma2^2*F11'*F11;
    %Fpinv = inv(E'*E + gamma^2*eye(n_modes))*E';   
    Fpinv = (E'*E + gamma^2*eye(n_modes))\E'; 
    x = (Fpinv*y)';
