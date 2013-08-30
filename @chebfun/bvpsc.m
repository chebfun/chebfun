function varargout = bvpsc(odefun, jacfun, bcmat, y0)
%BVPSC   Spectral collocation solution to a boundary-value problem.
%   Y = BVPSC(ODEFUN, JACFUN, BCMATS, Y0) solves a nonlinear boundary-value
%   problem using a spectral collocation method.
%
%   ODEFUN specifies the differential equation exactly as in BVP4C and BVP5C
%   with the 'Vectorized' option set. That is, if x is a row vector/CHEBFUN and
%   each row of y represents the values of one variable in the system,
%   ODEFUN(x, y) should return the derivative dy/dx, one row per variable.
%
%   JACFUN gives the Jacobian of the differential equation with respect to the
%   system variables. It has the same syntax as the 'FJacobian' option for BVP4C
%   and BVP5C. It is not vectorized, so JACFUN(x, y) should be an m-by-m matrix
%   in a system with m variables.
%
%   BCMAT gives three matrices describing the boundary condition in the form
%
%     BCMAT{1}*y(:,a) + BCMAT{2}*y(:,b) = BCMAT{3}.
%
%   The first two matrices are m-by-m and the third is an m-by-1 vector.
%   Despite the generality of the expression above, the boundary conditions
%   should be separated: for each i, row i of either BCMAT{1} or BCMAT{2}
%   must be zero.
%
%   Y0 should be a CHEBFUN that gives an initial guess to the solution. If the
%   residual of the equation with this guess is judged to be too large, BVP5C is
%   called for an improved guess, and it is initialized with Y0.
%
% See also BVP4C, BVP5C, BVPSET.

% [TODO]: Give an example?
% [TODO]: Does this method even work?

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

d = y0.domain;
x = chebfun(@(x) x, d);
if ( y0.isTransposed )
    y0 = y0.';
end

% If the initial residual is too big, use BVP5C() to initialize.
resid = diff(y0) - odefun(x.', y0.').';
if ( norm(resid, 'fro') > 1e-1*norm(y0, 'fro') )
    opt = bvpset('AbsTol', 1e-2, 'RelTol', 1e-2, 'Vectorized', 'on',...
        'fjacobian', jacfun);
    bcfun = @(ya, yb) bcmat{1}*ya + bcmat{2}*yb - bcmat{3};
    y = bvp5c(odefun, bcfun, y0, opt);
else
    y = y0;
end

% Set up some basic linear operator building blocks:
d = domain(d);
I = eye(d);
D = diff(d);
m = size(y, 2);
bigD = blkdiag(D, m);

state = warning('off', 'CHEBFUN:vectorwrap:shape');
ndy = zeros(8, 1);
% Newton iteration:
for newt = 1:8
    % Residual:
    resid = diff(y) - odefun(x', y.').';
    % Jacobian of the RHS odefun:
    J = jacobian(y);  
    bc = J.bc;
    J = (bigD - J) & bc;
    % Newton step:
    dy = J\resid;
    % Newton update:
    y = y - dy;
    % Stop if well converged, or stagnated.
    ndy(newt) = norm(dy, 'fro');
    if ( ndy(newt) < 1e-12 )
        break
    elseif ( newt > 3 && ndy(newt) > 0.1*ndy(newt-3) )
        warning('CHEBFUN:bvpsc', 'Newton iteration stagnated.')
        break
    end
end
warning(state)

% % Compress each variable to minimial representation.
% for j = 1:m
%     y(:,j) = chebfun(@(x) y(x,j), d);
% end

varargout{1} = y;

    function J = jacobian(y)
        % Jacobian due to the ODE. The Jacobian has diagonal blocks. The user
        % gives the entire J(x) for a single value of x. What we need is J_{i,j}
        % at all values of x to be encapsulated in a CHEBFUN. To minmize calls
        % to the jacfun, we do a construction of an arbitrary linear combination
        % of all elements, storing results in a 3D array. Then we go back and
        % convert each slice to a CHEBFUN.
        
        % randvec = rand(m^2,1);
        randvec = 1+.5*sin(1:m^2).';  % Arbitrary combination of variables
        A = [];
        function v = value(x)
            N = length(x);
            A = zeros(m, m, N);
            for l = 1:N
                A(:,:,l) = jacfun(x(l), feval(y, x(l)).');
            end
            v = (randvec' * reshape(A, [m^2, N])).';
        end
        
        chebfun(@value, d); % The random linear combination
        
        % Here each block is formed by CHEBFUN construction on the (i,j) element
        % of the Jacobian returned by the user's function. This means most of
        % the elements of the Jacobian are thrown away after each evaluation!
        J = chebop;
        for i = 1:m
            Jrow = chebop();
            for j = 1:m
                H = chebfun(A(i,j,:), d);
                Jrow = [Jrow, diag(H)];
            end
            J = [J ; Jrow];
        end
        
        % Attach the boundary conditions.
        Ka = bcmat{1};  
        Kb = bcmat{2};
        p = size(Ka,1);  
        kl = 1; 
        kr = 1;
        for k = 1:p
            atleft = any(Ka(k,:));
            atright = any(Kb(k,:));
            if ( ~xor(atleft, atright) )
                error('CHEBFUN:bvpsc:badBC',...
                    'Boundary conditions must be of separated type.')
            end
            bc = chebop;
            if ( atleft )
                for j = 1:m
                    bc = [bc, Ka(k,j)*I];
                end
                J.lbc(kl) = bc; 
                kl = kl + 1;
            else
                for j = 1:m
                    bc = [bc, Kb(k,j)*I];
                end
                J.rbc(kr) = bc;
                kr = kr  +1;
            end
        end
    end


end
