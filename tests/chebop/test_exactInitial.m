function pass = test_exactInitial(pref)
% Test that we are able to deal with a user passing the solution as the initial
% guess to a nonlinear BVP.

if ( nargin == 0 )
    pref = cheboppref;
end
tol = 1e-9;

try
    %% Problem set-up
    % Define the domain.
    dom = [0 10];
    
    % Assign the differential equation to a chebop on that domain.
    N = chebop(@(x,u) diff(u,2)+sin(u), dom);
    
    % Set up the rhs of the differential equation so that N(u) = rhs.
    rhs = 0;
    
    % Assign boundary conditions to the chebop.
    N.bc = @(x,u) [u(0)-2; u(10)-2];
    
    % Construct a linear chebfun on the domain,
    x = chebfun(@(x) x, dom);
    % and assign an initial guess to the chebop.
    N.init =  2.*cos(2.*pi.*x./10);
        
    % Option for discretization (either 'collocation' or 'ultraspherical').
    pref.discretization = @colloc2;
    
    % Solve to obtain a solution:
    u = solvebvp(N, rhs, pref);
    
    % Solve again, using the above solution as the initial guess
    N.init = u;
    u = solvebvp(N, rhs, pref);
    
    % Solve again, using a different discretization
    pref.discretization = @colloc1;
    u = solvebvp(N, rhs, pref);
     
    % Solve again, using a different discretization
    pref.discretization = @ultraS;
    u = solvebvp(N, rhs, pref);
    
    % Check that we are still happy with the solution
    err = norm(N(u));
    
    pass = err < tol; 
catch
    % Something went wrong, we've failed.
    pass = 0;
end