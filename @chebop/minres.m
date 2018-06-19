function [u,flag,relres,iter,resvec] = minres(L,f,tol,maxit,R1,R2,u0,varargin)
%MINRES    A Preconditioned Minimal Residual method for ODEs.
%   U = MINRES(L,F) attempts to solve the linear ODE L(U) = F on DOMAIN(L) with
%   boundary conditions. The CHEBOP L must be a self-adjoint, uniformly elliptic
%   second-order differential operator of the form
%          L(U) = (a(x)*U')' + c(x)U,
%   with dirichlet boundary conditions. The righthand side F should be a chebfun
%   on DOMAIN(L). By default the indefinite integral operator is employed as a
%   preconditioner.
%
%   U = MINRES(L,F,TOL) specifies the tolerance of the method. If TOL is [] then
%   MINRES uses the default in cheboppref.
%
%   U = MINRES(L,F,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then MINRES uses the default in cheboppref.
%
%   U = MINRES(L,F,TOL,MAXIT,R1,R2) solves the preconditioned linear ODE of
%   (R2*L*R1)(V) = R2*f, where R2 must be the adjoint of R1. R1 and R2 must be
%   function handles. If R1 = [], then the default preconditioner is employed.
%   Only the default preconditioner is currently supported.
%
%   U = MINRES(L,F,TOL,MAXIT,R1,R2,U0) specifies the initial guess. If U0 is []
%   then MINRES MINRES the default, the zero function on DOMAIN(L).
%
%   [U,FLAG] = MINRES(L,F,...) also returns a convergence FLAG:
%    0 MINRES converged to the desired tolerance TOL within MAXIT iterations
%    1 MINRES iterated MAXIT times but did not converge.
%    2 preconditioner R1 was an unbounded operator.
%    3 MINRES stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during MINRES became too
%      small or too large to continue computing.
%
%   [U,FLAG,RELRES] = MINRES(L,F,...) also returns the relative residual
%    NORM(R2(F)-R2(L(U)),2)/NORM(F,2). If FLAG is 0, then RELRES <= TOL.
%
%   [U,FLAG,RELRES,ITER] = MINRES(L,F,...) also returns the iteration number at
%   which U was computed: 0 <= ITER <= MAXIT.
%
%   [U,FLAG,RELRES,ITER,RESVEC] = MINRES(L,F,...) returns a vector of estimated
%   residual norms at each iteration including NORM(R2(F)-R2(L(Uk)),2).
%
% See also CHEBOP/PCG and CHEBOP/GMRES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Only continue if L is a linear chebop:
if ( ~all(islinear(L)) )
    error('CHEBFUN:CHEBOP:pcg:nonlinear', ...
        'MINRES supports only linear CHEBOP instances.');
end

% At the moment, we need a second-order differential equation:
LinearOp = linop( L );
if ( LinearOp.diffOrder ~=2 )
    error('CHEBFUN:CHEBOP:pcg:DiffOrder', ...
        'MINRES supports only second-order ODEs.');
end

if ( nargin < 2 )
    error(message('chebop:minres:NotEnoughInputs'));
end

% Norm of righthand side function:
n2f = norm(f, 2);
n = length( f );
dom = domain( f );

% Grab boundary conditions from chebop: 
if ( isa( L, 'chebop' ) )
    if isempty(L.lbc)
        left_bc = 0;
    else
        left_bc = L.lbcShow;
        
        % Current implementation requires Dirichlet boundary conditions: 
        if ( ~isa(left_bc, 'double' ) )
            error('CHEBFUN:CHEBOP:pcg:leftbc', ...
                   'Currently, we require Dirichlet boundary conditions. Please supply N.lbc = double.');
        end
    end
    
    if ( isempty(L.rbc) )
        right_bc = 0;
    else
        right_bc = L.rbcShow;
        % Current implementation requires Dirichlet boundary conditions: 
        if ( ~isa(right_bc, 'double' ) )
            error('CHEBFUN:CHEBOP:pcg:rightbc', ...
                   'Currently, we require Dirichlet boundary conditions. Please supply N.rbc = double.');
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

% Data mine L for faster "mat-vec" computations in the CG iteration:
x = chebfun( @(x) x );
c = L(x, 1+0*x);
a = -L(x, x.^2/2) - (-L(x, x) + c.*x).*x + c.*(x.^2/2);
L = @(v) -diff( a.*diff( v ) ) + c.*v;

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
if ( (nargin < 4) || isempty(maxit) )
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
else
    % The method only works with the default preconditioner:
    error(message('chebop:pcg:OnlyDefaultPreconditionerAllowed'))
end

% Projection operator
Pi = @(g) g - mean(g);

% Build modified differential operator:
x = chebfun( @(x) x, dom );
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

if ( nargin > 7 )
    error(message('chebop:minres:TooManyInputs'));
end

% Ensure that rhs is within the correct space, and if not then solve a modified
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
    warning('on','all');
    % g is now in the space Wn = {v in L_2: R1(v) in V_{n,0} }:
    g = Pi( R2f - R2( L(z) ) );
else
    % f is in the space Wn = {v in L_2: R1(v) in V_{n,0} }:
    g = PiR2f;
    z = 0*f;
end

% Set up for the method
flag = 1;
iter = 0;
umin = u;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolg = tol * normest(g);           % Relative tolerance
r = g - Tu;
normr = norm(r,2);                 % Norm of residual
normr_act = normr;

if ( normr <= tolg )               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2f;
    resvec = normr;
    % Undo preconditioner:
    u = R1( Pi( u ) );
    return
end

resvec = zeros(maxit+1,1);         % Preallocate vector for MINRES residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of minimum residual

vold = r;
v = vold;
beta1 = innerProduct(vold, v);
if ( beta1 <= 0 )
    flag = 5;
    relres = normr / n2b;
    resvec = resvec(1);
    u = R1( Pi( u ) );
    return
end
beta1 = sqrt(beta1);
snprod = beta1;
vv = v / beta1;
v = T(vv);
Amvv = v;
alpha = innerProduct(vv,v);
v = v - (alpha/beta1) * vold;

% Local reorthogonalization
numer = innerProduct(vv, v);
denom = innerProduct(vv, vv);
v = v - (numer/denom) * vv;
volder = vold;
vold = v;
betaold = beta1;
beta = innerProduct(v,v);
if ( beta < 0 )
    flag = 5;
    relres = normr / n2b;
    resvec = resvec(1);
    u = R1( Pi( u ) );
    return
end
iter = 1;
beta = sqrt(beta);
gammabar = alpha;
epsilon = 0;
deltabar = beta;
gamma = sqrt(gammabar^2 + beta^2);
mold = chebfun('0');
Amold = mold;
m = vv / gamma;
Am = Amvv / gamma;
cs = gammabar / gamma;
sn = beta / gamma;
u = u + snprod * cs * m;
snprod = snprod * sn;

normr = abs(snprod);
resvec(2,1) = normr;

% Check for convergence after first step.
if ( normr <= tolg )
    flag = 0;
    relres = normr / n2f;
    resvec = resvec(1:2);
    u = R1( Pi( u ) );
    return
end

stag = 0;                          % stagnation of the method
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)
for ii = 2:maxit
    
    vv = v * (1/beta);
    v = T(vv);
    Amolder = Amold;
    Amold = Am;
    Am = v;
    v = v - (beta / betaold) * volder;
    alpha = innerProduct(vv, v);
    v = v - (alpha / beta) * vold;
    volder = vold;
    vold = v;
    betaold = beta;
    beta = innerProduct(v, v);
    if (beta < 0)
        flag = 5;
        break
    end
    beta = sqrt(beta);
    delta = cs * deltabar + sn * alpha;
    molder = mold;
    mold = m;
    m = vv - delta * mold - epsilon * molder;
    Am = Am - delta * Amold - epsilon * Amolder;
    gammabar = sn * deltabar - cs * alpha;
    epsilon = sn * beta;
    deltabar = - cs * beta;
    gamma = sqrt(gammabar^2 + beta^2);
    m = m / gamma;
    Am = Am / gamma;
    cs = gammabar / gamma;
    sn = beta / gamma;
    % Check for stagnation of the method
    if (snprod*cs == 0) || (abs(snprod*cs)*norm(m) < eps*norm(u))
        % increment the number of consecutive iterates which are the same
        stag = stag + 1;
    else
        stag = 0;
    end
    u = u + (snprod * cs) * m;
    snprod = snprod * sn;
    normr = abs(snprod);
    
    resvec(ii+1,1) = normr;
    
    % Check for convergence
    if ( ( normr <= tolg ) || ( stag >= maxstagsteps ) || moresteps)
        % double check residual norm is less than tolerance
        r = g - T( u );
        normr_act = norm(r);
        resvec(ii+1, 1) = normr_act;
        if (normr_act <= tolg)
            flag = 0;
            iter = ii;
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning('tooSmallTolerance');
                end
                flag = 3;
                iter = ii;
                break;
            end
        end
    end
    
    if ( normr < normrmin )        % update minimal norm quantities
        normrmin = normr;
        umin = u;
        imin = ii;
    end
    
    if ( stag >= maxstagsteps )    % 3 iterates are the same
        flag = 3;
        break
    end
end                                % for ii = 1 : maxit
if ( isempty(ii) )
    ii = 1;
end

% Returned solution is first with minimal residual:
if ( flag == 0 )
    relres = normr_act / n2f;
else
    r_comp = g - T(u);
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

% truncate the zeros from resvec
if ( (flag <= 1) || (flag == 3) )
    resvec = resvec(1:ii+1);
else
    resvec = resvec(1:ii);
end

end