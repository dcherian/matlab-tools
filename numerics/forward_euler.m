% Implements explicit Forward Euler timestepping for parabolic equation
%       du/dt = d/dz(flux(u,z,t)) + source(u,z,t)
%       [] = forward_euler(fluxfn,sourcefn,sol)
%           fluxfn - function handle for flux term
%           sourcefn - function handle for source term
%           sol - structure with
%               - dt    - timestep
%               -  z    - grid
%               - nt    - number of timesteps
%               - uinit - initial guess

function [u] = forward_euler(fluxfn,sourcefn,sol)

    u = nan(length(uinit),nt); % u(z,t)
    try
        u(:,1) = sol.uinit;
    catch ME
        u(:,1) = sol.uinit';
    end
    
    dz = diff(z);

    for ii = 2:sol.nt
        flux = fluxfn(u(:,ii-1));
        source = sourcefn(u(:,ii-1));
        u(:,ii) = u(:,ii-1) ...
                    + sol.dt ./dz .* (diff(flux,1,1)) + source;
    end

    