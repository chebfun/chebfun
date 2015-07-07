function [u,lam,iter, retract] = newtonBVP(H,uinit, laminit, udot,lamdot, prefs)

% Find a tangent to the curve H(u,lambda)=0 at the given point, by solving
% a boundary-value problem.

d = domain(udot);
x = chebfun(@(x) x,d);
% Begin by creating a chebarray for the constraint which goes at the
% bottom
Jm = [udot'*diag(sum(d)) lamdot];

% Convert lam to chebconst, and store udot as u for when we iterate
lam = laminit;   % to get autodiff data
u = uinit;
% Compute Newton correction
accept = 0; iter = 0; retract = 0;
while ~accept
    % Evaluate residual, and linearise
    [S, res] = linearize(H, [u; lam], x);
    Scon = S.constraint;
    Scon.values = 0*Scon.values;
    
    % Create a chebarray
    L = linop([S; Jm]);
    L.constraint = Scon;
        
    
    % rhs is the residual, with zero at the bottom for the functional
    % condition
    rhs = [-res{1};0];
        
    dudlam = linsolve(L, rhs, prefs);
    du = dudlam{1};
    dlam = dudlam{2};
    
    % Should be close to solution, so take full Newton steps
    u = u+du;
    lam = lam+dlam;

    iter = iter + 1;
    if norm(du,2) < 1e-3
        accept = 1;
    elseif iter >=5 % Too many iterations
        retract = 1; return
    end

end
end