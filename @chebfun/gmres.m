function [u, flag, normres, Q] = gmres(varargin)
%GMRES   Iterative solution of chebfun operator equations.
%   U = GMRES(A, F) attempts to solve the operator equation A(U) = F, where F
%   and U are CHEBFUNs and A is a function handle defining a linear operator
%   on CHEBFUN.
%
%   U = GMRES(A, F, RESTART) chooses a restart parameter. Use Inf or [] for no
%   restarts. The default is Inf.
%
%   U = GMRES(A, F, RESTART, TOL) tries to find a solution for which the
%   residual norm is less than TOL relative to the norm of F. The default
%   value is 1e-10.
%
%   U = GMRES(A, F, RESTART, TOL, MAXIT) stops after MAXIT total iterations.
%   This defaults to 36.
%
%   U = GMRES(A, F, RESTART, TOL, MAXIT, M1INV, M2INV) applies the left and
%   right preconditioners M1INV and M2INV, defined as functions. They should
%   approximate the inverse of A. Note that this usage of preconditioners
%   differs from that in the built-in GMRES.
%
%   [U, FLAG] = GMRES(A,F,...) also returns a convergence FLAG:
%     0 GMRES converged to the desired tolerance TOL within MAXIT iterations.
%     1 GMRES iterated MAXIT times but did not converge.
%
%   [U, FLAG, NORMRES] = GMRES(A, F, ...) also returns a vector of the relative
%   residual norms for all iterations. Note the output ordering is not the same
%   as for built-in GMRES. This calling sequence will also print out updates on
%   the progress of the iteration.
%
% Example:
%   % To solve a simple Volterra integral equation:
%   f = chebfun('exp(-4*x.^2)', [-1 1]);
%   A = @(u) cumsum(u) + 20*u;
%   u = gmres(A, f, Inf, 1e-14);
%
% See also GMRES, CHEBOP/MLDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs and supply defaults:
defaults = {[], [], Inf, 1e-10, 36, [], []};
idx = nargin+1:length(defaults);
args = [varargin, defaults(idx)];
[L, f, m, tol, maxiter, M1inv, M2inv] = deal( args{:} );
showtrace = (nargout > 2);

if ( m == 0 )
    % No restarts:
    m = Inf;
end
% Avoid warning about Inf in for loop:
m = min(m, maxiter + 1);

% Initialise:
u = chebfun(0, domain(f));
normb = norm(f);
r = f;
if ( ~isempty(M1inv) )
    r = M1inv(r);
end
normres(1) = norm(r)/normb;

j = 1;                                             % total iterations
while ( (normres(j) > tol) && (j < maxiter) )      % outer iterations
    
    H = [];
    QTb = norm(r);
    Q = r/QTb ;                                    % Krylov basis
    
    for n = 1:m                                    % inner iterations
        
        % Next Krylov vector, with preconditioners.
        q = Q(:,n);
        if ( ~isempty(M2inv) )
            q = M2inv(q);
        end
        v = L(q);
        
        if ( ~isempty(M1inv) )
            v = M1inv(v);
        end
        
        % Modified Gram-Schmidt iteration.
        for k = 1:n
            H(k,n) = Q(:,k)' * v ;                   %#ok<AGROW>
            v = v - H(k,n)*Q(:,k);
        end
        H(n+1,n) = norm(v);                          %#ok<AGROW>
        
        % Use QR factorization to find the residual norm.
        % TODO: This could be made more efficient (worthwhile?).
        QTb(n+1,1) = 0;                              % by orthogonality
        j = j + 1;
        [P, R] = qr(H);
        normres(j) = abs( P(1,n+1)*QTb(1) ) / normb; %#ok<AGROW>
        
        % Done?
        if ( normres(j) < tol )
            showtrace = false;
            flag = 0;
            break
        end
        
        % Trace iteration--for debugging, not an official option.
        if ( showtrace && (rem(j - 1, 5) == 0) )
            fprintf('Iteration %2i: relative residual norm = %.3e\n', ...
                j - 1, normres(j));
        end
        
        % Give up?
        if ( j == maxiter + 1 )
            warning('CHEBFUN:CHEBFUN:gmres:maxiter', ...
                'Maximum number of iterations reached.')
            showtrace = false;
            break
        end
        
        % New basis vector:
        Q(:,n+1) = v / H(n+1,n);                    
        
        % Reorthogonalize (for research only--not an official option).
        if ( rem(n, Inf) == 0 )
            [Q, ignored] = qr(Q, 0);
        end
        
    end   % end inner iterations
    
    y = R(1:n,1:n)\(P(:,1:n)'*QTb);                % least squares soln
    u0 = Q(:,1:n)*y;                               % new part of solution
    if ( ~isempty(M2inv) )
        u0 = M2inv(u0);
    end
    u = u + u0;                                    % solution
    v = L(u0);
    if ( ~isempty(M1inv) )
        v = M1inv(v);
    end
    r = r - v;                                     % new residual
    if ( showtrace )
        fprintf('  (restart)\n')
    end
    
end   % end outer iterations

if ( j >= maxiter )
    flag = 1;
end
if ( nargout < 2 )
    fprintf('\n  Final relative residual: %.3e\n\n', normres(end))
end

end   % end main function

