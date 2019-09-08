function [u,flag,relres,iter,resvec] = gmres(L,f,restart,tol,maxit,R1,R2,u0,varargin)
%GMRES    A Preconditioned Generalized Minimal Residual method for ODEs.
%   U = GMRES(L,F) attempts to solve the linear ODE L(U) = F on DOMAIN(L) with
%   boundary conditions. The CHEBOP L must be a self-adjoint, uniformly elliptic
%   second-order differential operator of the form
%          L(U) = (a(x)*U')' b(x)*U' + c(x)U,
%   with dirichlet boundary conditions. The righthand side F should be a chebfun
%   on DOMAIN(L). By default the indefinite integral operator is employed as a
%   preconditioner. This uses the unrestarted method with the maximum number of
%   iteration set in chebfun.pref.
%
%   U = GMRES(L,F,RESTART) restarts the method every RESTART iterations. If
%   RESTART is N or [] then GMRES uses the unrestarted method as above.
%
%   U = GMRES(L,F,RESTART,TOL) specifies the tolerance of the method.  If TOL is
%   [] then GMRES uses the default, cheboppref().bvpTol.
%
%   U = GMRES(L,F,RESTART,TOL,MAXIT) specifies the maximum number of outer
%   iterations. Note: the total number of iterations is RESTART*MAXIT. If MAXIT
%   is [] then GMRES uses the default,  cheboppref().maxIter; If RESTART is []
%   then the total number of iterations is cheboppref().maxIter.
%
%   U = GMRES(L,F,RESTART,TOL,MAXIT,R1,R2) use preconditioner M or M=M1*M2 and
%   effectively solve the system R2(L(R1(U))) = R2(F) for U. Only the default
%   preconditioner is currently supported.
%
%   U = GMRES(L,F,RESTART,TOL,MAXIT,R1,R2,U0) specifies the initial guess. If U0
%   is [] then GMRES uses the default, the zero function on DOMAIN(L).
%
%   [U,FLAG] = GMRES(L,F,...) also returns a convergence FLAG:
%    0 GMRES converged to the desired tolerance TOL within MAXIT iterations.
%    1 GMRES iterated MAXIT times but did not converge.
%    2 preconditioner R1 was an unbounded operator.
%    3 GMRES stagnated (two consecutive iterates were the same).
%
%   [U,FLAG,RELRES] = GMRES(L,F,...) also returns the relative residual
%   NORM(R2(F)-R2(L(U)))/NORM(R2(F)). If FLAG is 0, then RELRES <= TOL.
%
%   [U,FLAG,RELRES,ITER] = GMRES(L,F,...) also returns both the outer and inner
%   iteration numbers at which X was computed: 0 <= ITER(1) <= MAXIT and 0 <=
%   ITER(2) <= RESTART.
%
%   [U,FLAG,RELRES,ITER,RESVEC] = GMRES(L,F,...) also returns a vector of the
%   residual norms at each inner iteration, including NORM(R2(F)-R2(L(U0)).
%
% See also CHEBOP/PCG and CHEBOP/MINRES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Only continue if L is a linear chebop:
if ( ~all(islinear(L)) )
    error('CHEBFUN:CHEBOP:pcg:nonlinear', ...
        'GMRES supports only linear CHEBOP instances.');
end

% At the moment, we need a second-order differential equation:
LinearOp = linop( L );
if ( LinearOp.diffOrder ~= 2 )
    error('CHEBFUN:CHEBOP:pcg:DiffOrder', ...
        'GMRES supports only second-order ODEs.');
end

if (nargin < 2)
    error(message('chebop:gmres:NotEnoughInputs'));
end

% Norm of righthand side function:
n2f = norm(f, 2);
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
                   'GMRES only supports Dirichlet boundary conditions. Please supply N.lbc = double.');
        end
    end
    
    if isempty(L.rbc)
        right_bc = 0;
    else
        right_bc = L.rbcShow;
        % Current implementation requires Dirichlet boundary conditions: 
        if ( ~isa(right_bc, 'double' ) )
            error('CHEBFUN:CHEBOP:pcg:rightbc', ...
                   'GMRES only supports Dirichlet boundary conditions. Please supply N.lbc = double.');
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
% Uses a non-divergence form
x = chebfun(@(x) x, dom);
c = L(x, 1+0*x); 
bminusa = L(x, x) - c.*x;
a = -L(x, x.^2/2) + bminusa.*x + c.*x.^2/2;
b = bminusa + diff(a);
L = @(v) -diff( a.*diff( v ) ) + b.*diff( v ) +  c.*v;

% Assign default values to unspecified parameters
if ( (nargin < 3) || isempty(restart)  )
    restarted = false;
else
    restarted = true;
end

if ( (nargin < 4) || isempty(tol) )
    tol = cheboppref().bvpTol;
end

warned = 0;
if ( tol < eps )
    warning(message('chebop:gmres:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif ( tol >= 1 )
    warning(message('chebop:gmres:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end

if (nargin < 5) || isempty(maxit)
    maxit = cheboppref().maxIter;
end

if ( restarted )
    outer = maxit;
    inner = restart;
else
    outer = 1;
    inner = maxit;
end

if ( (nargin < 6) || isempty(R1) )
    % Preconditioner:
    R1 = @(u) cumsum(u);
else
   % The method only works with the default preconditioner:
   error(message('chebop:gmres:OnlyDefaultPreconditionerAllowed'))
end

if ( (nargin < 7) || isempty(R2) )
    % Adjoint operator to R1:
    R2 = @(u) sum(u) - cumsum(u);
else
   % The method only works with the default preconditioner:
   error(message('chebop:gmres:OnlyDefaultPreconditionerAllowed'))
end

% Projection operator
Pi = @(g) g - mean(g);

% Build modified differential operator:
T = @(v) Pi( R2( L( R1( v ) ) ) );

if ((nargin >= 8) && ~isempty(u0))
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

if (nargin > 8) 
    error(message('chebop:gmres:TooManyInputs'));
end

% Ensure that rhs is within the correct space, and if not then solve a
% modified problem. 
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
umin = u;                        % Iterate which has minimal residual so far
imin = 0;                        % "Outer" iteration at which xmin was computed
jmin = 0;                        % "Inner" iteration at which xmin was computed
tolg = tol * norm(g);                % Relative tolerance
stag = 0;
moresteps = 0;
maxmsteps = 5;
maxstagsteps = 3;
minupdated = 0;

r = g - Tu;
normr = norm(r,2);               % Norm of residual
normr_act = normr;

if (normr <= tolg)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2f;
    resvec = normr;
    % Undo preconditioner:
    u = R1( Pi( u ) ) + z;
    return
end

normr = norm(r);                 % norm of the preconditioned residual
n2g = norm(g);                   % norm of the preconditioned rhs
tolg = tol * n2g;
if (normr <= tolg)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2g;
    iter = [0 0];
    resvec = n2g;
    return
end

resvec = zeros(inner*outer+1,1); % Preallocate vector for norm of residuals
resvec(1) = normr;               % resvec(1) = norm(b-A*x0)
normrmin = normr;                % Norm of residual from xmin
R = zeros(inner,inner);

for outiter = 1:outer

    H = [];
    QTb = norm(r);
    Q = r/QTb ;                  % Krylov basis
    
    for initer = 1:inner
        
        % Next Krylov vector, with preconditioners.
        q = Q(:,initer);
        v = T(q);
        
        % Modified Gram-Schmidt iteration.
        for k = 1:initer
            H(k,initer) = Q(:,k)' * v ;                   
            v = v - H(k,initer)*Q(:,k);
        end
        H(initer+1,initer) = norm(v);                          
        
        % Use QR factorization to find the residual norm.
        % TODO: This could be made more efficient (worthwhile?).
        QTb(initer+1,1) = 0;      % by orthogonality
        [P, R] = qr(H);
 
        % New basis vector:
        Q(:,initer+1) = v / H(initer+1,initer);                    

        normr = abs( P(1,initer+1)*QTb(1) ); 
        resvec((outiter-1)*inner+initer+1) = normr;
        normr_act = normr;
        
        if ( normr <= tolg || stag >= maxstagsteps || moresteps )
            % Evalute this outer loop's solution
            y = R(1:initer,1:initer)\(P(:,1:initer)'*QTb);  % least squares soln
            additive = Q(:,1:initer)*y;                     % new part of solution
            if norm(additive) < eps*norm(u)
                stag = stag + 1;
            else
                stag = 0;
            end
            % Current total solution
            um = u + additive;

            % Evaluate current residual
            r = g - T(um);
            normr_act = norm(r);
            resvec((outiter-1)*inner+initer+1) = normr_act;
            
            if ( normr_act <= normrmin )
                normrmin = normr_act;
                imin = outiter;
                jmin = initer;
                umin = um;
                minupdated = 1;
            end
            
            if ( normr_act <= tolg )
                u = um;
                flag = 0;
                iter = [outiter, initer];
                break
            else
                if ( stag >= maxstagsteps && moresteps == 0 )
                    stag = 0;
                end
                moresteps = moresteps + 1;
                if ( moresteps >= maxmsteps )
                    if ( ~warned )
                        warning(message('chebop:gmres:tooSmallTolerance'));
                    end
                    flag = 3;
                    iter = [outiter, initer];
                    break;
                end
            end
        end
        
        if ( normr_act <= normrmin )
            normrmin = normr_act;
            imin = outiter;
            jmin = initer;
            minupdated = 1;
        end
        
        if ( stag >= maxstagsteps )
            flag = 3;
            break;
        end
    end         % ends inner loop

    if ( flag ~= 0 )
        if minupdated
            idx = jmin;
        else
            idx = initer;
        end
        y = R(1:idx,1:idx)\(P(:,1:idx)'*QTb); % least squares soln
        additive = Q(:,1:idx)*y;   
        u = u + additive;
        umin = u;
        r = g - T(u);
        normr_act = norm(r);
    end
    
    if ( normr_act <= normrmin )
        umin = u;
        normrmin = normr_act;
        imin = outiter;
        jmin = initer;
    end
    
    if ( flag == 3 )
        break
    end
    if ( normr_act <= tolg )
        flag = 0;
        iter = [outiter, initer];
        break
    end
    minupdated = 0;
end         % ends outer loop

% returned solution is that with minimum residual
if ( flag == 0 )
    relres = normr_act / n2g;
else
    u = umin;
    iter = [imin jmin];
    relres = normr_act / n2g;
end

%Undo preconditioner, read BC subsolution
u = R1(u) + z;

resvec = resvec(1:(outiter-1)*inner+initer+1);
if ( flag == 2 && initer ~= 0 )
    resvec(end) = [];
end

end
