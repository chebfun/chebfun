function varargout = gmres(A, varargin)
%GMRES   Iterative solution of a linear system. 
%   U = GMRES(A,F) solves the system A*U = F for CHEBFUN U and F and linear 
%   CHEBOP A. If A is not linear, an error is returned.
%
%   More calling options are available; see chebfun/gmres for details.
%
% Example:
%   % To solve a simple Volterra integral equation:
%   d = [-1,1];
%   f = chebfun('exp(-4*x.^2)',d);
%   A = chebop(@(u) cumsum(u) + 20*u, d);
%   u = gmres(A,f,Inf,1e-14);
%
% See also CHEBFUN/GMRES, GMRES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~all(islinear(A)) )
    error('CHEBFUN:CHEBOP:gmres:nonlinear', ...
        'GMRES supports only linear CHEBOP instances.');
end

% If a second-order differential operator, then call operator GMRES. 
% Otherwise, go to chebfun/gmres:
LinearOp = linop( A );
if ( LinearOp.diffOrder ==2 )
    [varargout{1:nargout}] = OperatorGMRES(A,varargin{:});
else
    op = @(u) A*u;
    [varargout{1:nargout}] = gmres(op, varargin{:});
end

end   

function [u,flag,relres,iter,resvec] = OperatorGMRES(L,f,restart,tol,maxit,R1,R2,u0,varargin)
%GMRES   Generalized Minimum Residual Method.
%   X = GMRES(A,B) attempts to solve the system of linear equations A*X = B
%   for X.  The N-by-N coefficient matrix A must be square and the right
%   hand side column vector B must have length N. This uses the unrestarted
%   method with MIN(N,10) total iterations.
%
%   X = GMRES(AFUN,B) accepts a function handle AFUN instead of the matrix
%   A. AFUN(X) accepts a vector input X and returns the matrix-vector
%   product A*X. In all of the following syntaxes, you can replace A by
%   AFUN.
%
%   X = GMRES(A,B,RESTART) restarts the method every RESTART iterations.
%   If RESTART is N or [] then GMRES uses the unrestarted method as above.
%
%   X = GMRES(A,B,RESTART,TOL) specifies the tolerance of the method.  If
%   TOL is [] then GMRES uses the default, 1e-6.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT) specifies the maximum number of outer
%   iterations. Note: the total number of iterations is RESTART*MAXIT. If
%   MAXIT is [] then GMRES uses the default, MIN(N/RESTART,10). If RESTART
%   is N or [] then the total number of iterations is MAXIT.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M) and
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M1,M2) use preconditioner M or M=M1*M2
%   and effectively solve the system inv(M)*A*X = inv(M)*B for X. If M is
%   [] then a preconditioner is not applied.  M may be a function handle
%   returning M\X.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M1,M2,X0) specifies the first initial
%   guess. If X0 is [] then GMRES uses the default, an all zero vector.
%
%   [X,FLAG] = GMRES(A,B,...) also returns a convergence FLAG:
%    0 GMRES converged to the desired tolerance TOL within MAXIT iterations.
%    1 GMRES iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 GMRES stagnated (two consecutive iterates were the same).
%
%   [X,FLAG,RELRES] = GMRES(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL. Note with
%   preconditioners M1,M2, the residual is NORM(M2\(M1\(B-A*X))).
%
%   [X,FLAG,RELRES,ITER] = GMRES(A,B,...) also returns both the outer and
%   inner iteration numbers at which X was computed: 0 <= ITER(1) <= MAXIT
%   and 0 <= ITER(2) <= RESTART.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = GMRES(A,B,...) also returns a vector of
%   the residual norms at each inner iteration, including NORM(B-A*X0).
%   Note with preconditioners M1,M2, the residual is NORM(M2\(M1\(B-A*X))).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      x = gmres(A,b,10,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = afun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   and this preconditioner backsolve function
%      %------------------------------------------%
%      function y = mfun(r,n)
%      y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%      %------------------------------------------%
%   as inputs to GMRES:
%      x1 = gmres(@(x)afun(x,n),b,10,tol,maxit,@(x)mfun(x,n));
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, LSQR, MINRES, PCG, QMR, SYMMLQ,
%   TFQMR, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2015 The MathWorks, Inc.

if (nargin < 2)
    error(message('chebop:gmres:NotEnoughInputs'));
end

% Norm of righthand side function:
n2f = norm(f, 2);
dom = domain( f );

if ( isa( L, 'chebop' ) )
    if isempty(L.lbc)
        left_bc = 0;
    else
        left_bc = L.lbcShow;
    end
    if isempty(L.rbc)
        right_bc = 0;
    else
        right_bc = L.rbcShow;
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
x = chebfun( @(x) x );
c = L(x, 1+0*x); 
bminusa = L(x, x) - c.*x;
a = -L(x, x.^2/2) + bminusa.*x + c.*x.^2/2;
b = bminusa + diff(a);
L = @(v) -diff( a.*diff( v ) ) + b.*diff( v ) +  c.*v;

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(restart) 
    restarted = false;
else
    restarted = true;
end

if (nargin < 4) || isempty(tol)
    tol = cheboppref().bvpTol;
end

warned = 0;
if tol < eps
    warning(message('chebop:gmres:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('chebop:gmres:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end

if (nargin < 5) || isempty(maxit)
    maxit = cheboppref().maxIter;
end

if restarted
    outer = maxit;
    inner = restart;
else
    outer = 1;
    inner = maxit;
end


if ((nargin < 6) || isempty(R1))
    % Preconditioner:
    R1 = @(u) cumsum(u);
else
   % The method only works with the default preconditioner:
   error(message('chebop:gmres:OnlyDefaultPreconditionerAllowed'))
end

if ((nargin < 7) || isempty(R2))
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
evalxm = 0;
stag = 0;
moresteps = 0;
maxmsteps = 5;
maxstagsteps = 3;
minupdated = 0;

r = g - Tu;
normr = norm(r,2);                   % Norm of residual
normr_act = normr;

if (normr <= tolg)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2f;
    resvec = normr;
    % Undo preconditioner:
    u = R1( Pi( u ) ) + z;
    return
end


normr = norm(r);                 % norm of the preconditioned residual
n2g = norm(g);         % norm of the preconditioned rhs
tolg = tol * n2g;
if (normr <= tolg)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2g;
    iter = [0 0];
    resvec = n2g;
    return
end

resvec = zeros(inner*outer+1,1);  % Preallocate vector for norm of residuals
resvec(1) = normr;                % resvec(1) = norm(b-A*x0)
normrmin = normr;                 % Norm of residual from xmin


R = zeros(inner,inner);

for outiter = 1 : outer
    
    
    H = [];
    QTb = norm(r);
    Q = r/QTb ;                                    % Krylov basis
    
    for initer = 1 : inner
        
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
        QTb(initer+1,1) = 0;                              % by orthogonality
        [P, R] = qr(H);
 
        % New basis vector:
        Q(:,initer+1) = v / H(initer+1,initer);                    

        normr = abs( P(1,initer+1)*QTb(1) ); 
        resvec((outiter-1)*inner+initer+1) = normr;
        normr_act = normr;
        
        if (normr <= tolg || stag >= maxstagsteps || moresteps)
            % Evalute this outer loop's solution
            y = R(1:initer,1:initer)\(P(:,1:initer)'*QTb);                % least squares soln
            additive = Q(:,1:initer)*y;                               % new part of solution
            if norm(additive) < eps*norm(u)
                stag = stag + 1;
            else
                stag = 0;
            end
            % Current total solution
            um = u + additive;

            %Evaluate current residual
            r = g - T(um);
            normr_act = norm(r);
            resvec((outiter-1)*inner+initer+1) = normr_act;
            
            if normr_act <= normrmin
                normrmin = normr_act;
                imin = outiter;
                jmin = initer;
                umin = um;
                minupdated = 1;
            end
            
            if normr_act <= tolg
                u = um;
                flag = 0;
                iter = [outiter, initer];
                break
            else
                if stag >= maxstagsteps && moresteps == 0
                    stag = 0;
                end
                moresteps = moresteps + 1;
                if moresteps >= maxmsteps
                    if ~warned
                        warning(message('chebop:gmres:tooSmallTolerance'));
                    end
                    flag = 3;
                    iter = [outiter, initer];
                    break;
                end
            end
        end
        
        if normr_act <= normrmin
            normrmin = normr_act;
            imin = outiter;
            jmin = initer;
            minupdated = 1;
        end
        
        if stag >= maxstagsteps
            flag = 3;
            break;
        end
    end         % ends inner loop
    
    evalxm = 0;
    
    if flag ~= 0
        if minupdated
            idx = jmin;
        else
            idx = initer;
        end
        y = R(1:idx,1:idx)\(P(:,1:idx)'*QTb);                % least squares soln
        additive = Q(:,1:idx)*y;   
        u = u + additive;
        umin = u;
        r = g - T(u);
        normr_act = norm(r);
    end
    
    if normr_act <= normrmin
        umin = u;
        normrmin = normr_act;
        imin = outiter;
        jmin = initer;
    end
    
    if flag == 3
        break;
    end
    if normr_act <= tolg
        flag = 0;
        iter = [outiter, initer];
        break;
    end
    minupdated = 0;
end         % ends outer loop

% returned solution is that with minimum residual
if flag == 0
    relres = normr_act / n2g;
else
    u = umin;
    iter = [imin jmin];
    relres = normr_act / n2g;
end

%Undo preconditioner, readd BC subsolution
u = R1(u) + z;

resvec = resvec(1:(outiter-1)*inner+initer+1);
if flag == 2 && initer ~= 0
    resvec(end) = [];
end

end
