function [u,flag,relres,iter,resvec] = pcg(L,f,tol,maxit,R1,R2,u0,varargin)
%PCG    Preconditioned Conjugate Gradients Method for ODEs.
%   U = PCG(L,F) attempts to solve the linear ODE L(U) = F on DOMAIN(L) with
%   boundary conditions. The CHEBOP L must be a self-adjoint, uniformly elliptic
%   second-order differential operator of the form
%          L(U) = (a(x)*U')' + c(x)U,  c(x)>=0,
%   with dirichlet boundary conditions. The righthand side F should be a chebfun
%   on DOMAIN(L). By default the indefinite integral operator is employed as a
%   preconditioner.
%
%   U = PCG(L,F,TOL) specifies the tolerance of the method. If TOL is [] then
%   PCG uses the default in cheboppref.
%
%   U = pcg(L,F,TOL,MAXIT) specifies the maximum number of iterations. If MAXIT
%   is [] then pcg uses the default in cheboppref.
%
%   U = pcg(L,F,TOL,MAXIT,R1,R2) solves the preconditioned linear ODE of
%   (R2*L*R1)(V) = R2*f, where R2 must be the adjoint of R1. R1 and R2 must be
%   function handles. If R1 = [], then the default preconditioner is employed.
%   Currently only the default preconditioner is supported.
%
%   U = PCG(L,F,TOL,MAXIT,R1,R2,U0) specifies the initial guess. If U0 is []
%   then pcg uses the default, the zero function on DOMAIN(L).
%
%   [U,FLAG] = PCG(L,F,...) also returns a convergence FLAG:
%    0 pcg converged to the desired tolerance TOL within MAXIT iterations
%    1 pcg iterated MAXIT times but did not converge.
%    2 preconditioner R1 was an unbounded operator.
%    3 pcg stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during pcg became too
%      small or too large to continue computing.
%
%   [U,FLAG,RELRES] = PCG(L,F,...) also returns the relative residual
%   NORM(R2(F)-R2(L(U)), 2)/NORM(R2(F), 2). If FLAG is 0, then RELRES <= TOL.
%
%   [U,FLAG,RELRES,ITER] = PCG(L,F,...) also returns the iteration number at
%   which U was computed: 0 <= ITER <= MAXIT.
%
%   [U,FLAG,RELRES,ITER,RESVEC] = PCG(L,F,...) returns a vector of the estimated
%   residual norms at each iteration including NORM(R2(F)-R2(L(U0)), 2).
%
% See also CHEBOP/MINRES and CHEBOP/GMRES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Only continue if L is a linear chebop:
if ( ~all(islinear(L)) )
    error('CHEBFUN:CHEBOP:pcg:nonlinear', ...
        'PCG supports only linear CHEBOP instances.');
end

% At the moment, we need a second-order differential equation:
LinearOp = linop( L );
if ( LinearOp.diffOrder ~=2 )
    error('CHEBFUN:CHEBOP:pcg:DiffOrder', ....
        'PCG supports only second-order ODEs.');
end

if ( nargin < 2 )
    error(message('chebop:pcg:NotEnoughInputs'));
end

% Norm of righthand side function:
n2f = norm(f, 2);
dom = domain( f );

% Grab boundary conditions from chebop: 
if ( isa( L, 'chebop' ) )
    if ( isempty(L.lbc) )
        left_bc = 0;
    else
        left_bc = L.lbcShow;
        
        % Current implementation requires Dirichlet boundary conditions: 
        if ( ~isa(left_bc, 'double' ) )
            error('CHEBFUN:CHEBOP:pcg:leftbc', ...
                   'PCG only supports Dirichlet boundary conditions. Please supply N.lbc = double.');
        end
    end
    
    if ( isempty(L.rbc) )
        right_bc = 0;
    else
        right_bc = L.rbcShow;
        % Current implementation requires Dirichlet boundary conditions: 
        if ( ~isa(right_bc, 'double' ) )
            error('CHEBFUN:CHEBOP:pcg:rightbc', ...
                   'PCG only supports Dirichlet boundary conditions. Please supply N.rbc = double.');
        end
    end
    
    L = L.op;
elseif ( isa( L, 'function_handle' ) )
    left_bc = 0;
    right_bc = 0;
else
    error(message('chebop:pcg:DiffOperatorIllDefined'));
end

% Ensure that L = @(x,u) ...
if ( nargin(L) == 1 )
    L = @(x,u) L(u);
elseif ( nargin(L) ~=2 )
    error(message('chebop:pcg:DiffOpNargin'));
end

% Assign default values to unspecified parameters
if ( (nargin < 3) || isempty(tol) )
    tol = cheboppref().bvpTol;
end
warned = 0;
if ( tol <= eps )
    warning(message('chebop:pcg:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif ( tol >= 1 )
    warning(message('chebop:pcg:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 4) || isempty(maxit)
    maxit = cheboppref().maxIter;
end

if ( (nargin < 5) || isempty(R1) )
    % Preconditioner:
    R1 = @(u) cumsum(u);
else
    % The method only works with the default preconditioner:
    error(message('chebop:pcg:OnlyDefaultPreconditionerAllowed'))
end

if ( (nargin < 6) || isempty(R2) )
    % Adjoint operator to R1:
    R2 = @(u) sum(u) - cumsum(u);
end

% Projection operator
Pi = @(g) g - mean(g);

% Data mine L for faster "operator-function" products in the CG method:
x = chebfun( @(x) x, dom);
c = L(x, 1+0*x);
if ( min(c) < -tol )
    warning('Chebfun:chebop:pcg:Eigenvalues', ...
        'Does the differential operator have nonnegative eigenvalues?');
end
a = -L(x, x.^2/2) - (-L(x, x) + c.*x).*x + c.*(x.^2/2);
if ( min(a) < -tol )
    warning('Chebfun:chebop:pcg:elliptic', ...
        'Is the differential operator uniformly elliptic?');
end
L = @(v) -diff( a.*diff( v ) ) + c.*v;

% Build modified differential operator:
T = @(v) Pi( R2( L( R1( v ) ) ) );

if ( (nargin >= 7) && ~isempty(u0) )
    if ( ~domainCheck(f, u0) )
        error(message('chebop:pcg:WrongInitGuessDomain'));
    else
        u = u0;
    end
    Tu = T(u);
else
    % Otherwise, initial guess is the zero function:
    u = 0*f;
    Tu = u;
end

% Ensure that RHS is within the correct space, and if not then solve a modified
% problem.
R2f = R2( f );
PiR2f = Pi( R2f );
if ( norm( R2f - PiR2f ) > tol || norm(left_bc)> tol || norm(right_bc) > tol )
    % f is NOT in the space Wn = {v in L_2: R1(v) in V_{n,0} }:
    basis = x.^(0:4);
    A = zeros(4, size(basis,2));
    for jj = 1:size(basis,2)
        A(1:2,jj) = feval( R1( R2( L(basis(:,jj)) ) ), dom([1,length(dom)])');
        A(3:4,jj) = feval(basis(:,jj),dom([1,length(dom)])');
    end
    b = [ feval( R1( R2f ), dom([1,length(dom)])') ; left_bc ; right_bc];
    warning('off','all')
    z = basis*(A\b);
    warning('on','all')
    % g is now in the space Wn = {v in L_2: R1(v) in V_{n,0} }:
    g = Pi( R2f - R2( L(z) ) );
else
    % f is in the space Wn = {v in L_2: R1(v) in V_{n,0} }:
    g = PiR2f;
    z = 0*f;
end

% Set up for the method
flag = 1;
umin = u;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolf = tol*norm(g,2);              % Relative tolerance.
r = g - Tu;
p = r;
normr = norm(r,2);                 % Norm of residual
normr_act = normr;

% Initial guess is a good enough solution:
if ( normr <= tolf )
    flag = 0;
    relres = normr / n2f;
    iter = 0;
    resvec = normr;
    % Undo preconditioner:
    u = R1( Pi( u ) );
    return
end

resvec = zeros(maxit+1,1);   % Preallocate vector for norm of residuals
resvec(1,:) = normr;         % resvec(1) = norm(R2(f)-T(u0),2)
normrmin = normr;            % Norm of minimum residual
stag = 0;                    % stagnation of the method
moresteps = 0;
maxmsteps = 5;
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)
rho = innerProduct(r, r);
for ii = 1 : maxit
    
    % Analogue of mat-vec:
    Lp = T( p );
    
    % The preconditioned CG method for ODEs:
    alpha = rho./innerProduct(p, Lp);
    u = u + alpha*p;
    r = r - alpha*Lp;
    rho_new = innerProduct(r, r);
    beta = rho_new / rho;
    p = r + beta*p;
    rho = rho_new;
    
    % Check everything worked out OK:
    if ((rho == 0) || isinf(rho))
        flag = 4;   % Scalar quantity too small or too large
        normr_act = sqrt(rho);
        break
    end
    
    if ( isinf(alpha) )
        flag = 4;   % Scalar quantity too small or too large
        break
    end
    
    if ( (beta == 0) || isinf(beta) )
        flag = 4;   % Scalar quantity too small or too large
        break
    end
    
    % Check for stagnation of the method
    if ( norm(p)*abs(alpha) < eps*norm(u) )
        stag = stag + 1;
    else
        stag = 0;
    end
    
    % Store residual information:
    normr = sqrt(rho);
    normr_act = normr;
    resvec(ii+1, 1) = normr;
    
    % Check for convergence
    if ( ( normr <= tolf ) || ( stag >= maxstagsteps ) || moresteps )
        r = g - T( u );
        normr_act = norm(r);
        resvec(ii+1, 1) = normr_act;
        if ( normr_act <= tolf )
            flag = 0;     % successfully converged
            iter = ii;
            break
        else
            if ( ( stag >= maxstagsteps ) && ( moresteps == 0 ) )
                stag = 0;
            end
            moresteps = moresteps + 1;
            if ( moresteps >= maxmsteps )
                if ( ~warned )
                    warning('The tolerance is probably too small.');
                end
                flag = 3;
                iter = ii;
                break
            end
        end
    end
    
    % Update minimal norm quantities:
    if ( normr_act < normrmin )
        normrmin = normr_act;
        umin = u;
        imin = ii;
    end
    if ( stag >= maxstagsteps )
        flag = 3;     % PCG convergence has stagnated.
        break
    end
end % for ii = 1:maxit

% Returned solution is first with minimal residual:
if ( flag == 0 )
    relres = normr_act / n2f;
else
    r_comp = g - T(umin);
    if ( norm(r_comp) <= normr_act )
        u = umin;
        iter = imin;
        relres = norm(r_comp) / n2f;
    else
        iter = ii;
        relres = normr_act / n2f;
    end
end

% Undo preconditioner:
u = R1( u ) + z;

% Truncate zeros from resvec
if ( (flag <= 1) || (flag == 3) )
    resvec = resvec(1:ii+1,:);
else
    resvec = resvec(1:ii, :);
end

end