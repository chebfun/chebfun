function varargout = remez(f, varargin)
%REMEZ   Best polynomial or rational approximation for real valued chebfuns.
%   P = REMEZ(F, M) computes the minimax polynomial approximation of degree M
%   to the real CHEBFUN F using the Remez algorithm.
%
%   [P, Q] = REMEZ(F, M, N) computes the minimax rational approximation P/Q
%   of type (M, N).
%
%   [P, Q, R_HANDLE] = REMEZ(F, M, N) additionally returns a function handle
%   for evaluating P/Q.
%
%   [...] = REMEZ(..., 'tol', TOL) uses the value TOL as the termination
%   tolerance on the relative equioscillation error.
%
%   [...] = REMEZ(..., 'display', 'iter') displays output at each iteration.
%
%   [...] = REMEZ(..., 'maxiter', MAXITER) sets the maximum number of allowable
%   iterations to MAXITER.
%
%   [...] = REMEZ(..., 'plotfcns', 'error') plots the error after each iteration
%   while the algorithm executes.
%
%   [P, ERR] = REMEZ(...) and [P, Q, R_HANDLE, ERR] = REMEZ(...) returns the
%   maximum error ERR.
%
%   [P, ERR, STATUS] = REMEZ(...) and [P, Q, R_HANDLE, ERR, STATUS] = REMEZ(...)
%   return a structure array STATUS with the following fields:
%      STATUS.DELTA  - Obtained tolerance.
%      STATUS.ITER   - Number of iterations performed.
%      STATUS.DIFFX  - Maximum correction in last trial reference.
%      STATUS.XK     - Last trial reference on which the error equioscillates.
%
%   This code is quite reliable for polynomial approximations but may sometimes
%   have difficulties in the rational case.
%
% References:
%
%   [1] S. Filip and Y. Nakatsukasa, manuscript in preparation.
%
%   [2] R. Pachon and L. N. Trefethen, "Barycentric-Remez algorithms for best
%   polynomial approximation in the chebfun system", BIT Numerical Mathematics,
%   49:721-742, 2009.
%
%   [3] R. Pachon, "Algorithms for Polynomial and Rational Approximation".
%   D. Phil. Thesis, University of Oxford, 2010 (Chapter 6).
%
% See also CF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:remez:real', ...
        'REMEZ only supports real valued functions.');
end

if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:remez:quasi', ...
        'REMEZ does not currently support quasimatrices.');
end

if ( issing(f) )
    error('CHEBFUN:CHEBFUN:remez:singularFunction', ...
        'REMEZ does not currently support functions with singularities.');
end

% Parse the inputs.
[m, n, N, rationalMode, symFlag, xk, opts] = parseInputs(f, varargin{:});

% If m=-1, this means f=odd and input (m,n)=(0,n); return constant 0. 
if ( m == -1 )
    q = chebfun(1, f.domain([1,end]));
    p = chebfun(0, f.domain([1,end]));
    varargout = {p, q, @(x) feval(p, x)./feval(q, x), norm(f,'inf'), []};    
    return
end

if(isempty(xk)) % no initial reference is given by the user
    % Several initialization attempts are made
    if(n == 0) % polynomial case
        xk = chebpts(N + 2, f.domain([1, end]));
        [p,err,status] = remezKernel(f,m, n, N, rationalMode, xk, opts, 1);
        q = chebfun('1');
        varargout = {p,q,p,err,status};
    else % rational case
        % A first attempt is using CF as an initial guess
        disp('Trying CF-based initialization...');
        try
        xk = cfInit(f, m, n);
        [p,q,rh,err,status] = remezKernel(f,m, n, N, rationalMode, xk, opts, 1);
        varargout = {p,q,rh,err,status};
        catch
        disp('CF-based initialization failed, turning to AAA-Lawson')
        status.success = 0; % CF didn't work
        end
        
        % If CF doesn't give a satisfactory answer, we try AAA-Lawson
        if(status.success == 0)
            disp('Trying AAA-Lawson-based initialization...');
            xk = AAALawsonInit(f, m, n);
            [p,q,rh,err,status] = remezKernel(f, m, n, N, rationalMode, xk, opts, 1);
            varargout = {p,q,rh,err,status};
        end

        % A final attempt using cumulative distribution functions
        if(status.success == 0)
            xk = cdfInit(f,m,n,symFlag,opts,1);
            [p,q,rh,err,status] = remezKernel(f,m, n, N, rationalMode, xk, opts, 1);
            varargout = {p,q,rh,err,status};
        end
        
        if(status.success == 0)
            xk = cdfInit(f,m,n,symFlag,opts,2);
            [p,q,rh,err,status] = remezKernel(f,m, n, N, rationalMode, xk, opts, 1);
            varargout = {p,q,rh,err,status};
        end
    end
else  % the user has also given a starting reference
    if(n == 0)
        [p,err,status] = remezKernel(f,m, n, N, rationalMode, xk, opts, 1);
        varargout = {p,err,status};
    else
        [p,q,rh,err,status] = remezKernel(f,m, n, N, rationalMode, xk, opts, 1);
        varargout = {p,q,rh,err,status};
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implementing the different initialization strategies

% CF-based initialization
function xk = cfInit(f, m, n)
    if ( numel(f.funs) == 1 )
        [p, q] = cf(f, m, n);
    else
        [p, q] = cf(f, m, n, 50*(m+n) );
    end
    pqh = @(x) feval(p, x)./feval(q, x);
    [xk, ~, ~, flag] = exchange([], 0, 2, f, p, pqh, m+n+2, n);

    % If the above procedure failed to produce a reference
    % with enough oscillation points, use polynomial Remez.
    if ( flag == 0 )
        [~,~,~,~,status] = remez(f,m+n); xk = status.xk;
    end
end

% Now turn to initialization via AAA-Lawson, this is more expensive than CF
% but less so than CDF (which follows if this fails). 
% 
function xk = AAALawsonInit(f,m,n) % AAA-Lawson initialization for functions with breakpoints
    NN = max(10*max(m,n),round(1e5/max(m,n))); 
    dom = domain(f);
    Z = linspace(dom(1),dom(end),NN);   
    F = feval(f,Z);
    [r,~,~,~,xk] = aaamn_lawson(F,Z,m,n);    % 1st AAA-Lawson
    xk = findReference(f,r,m,n,xk);
    
    % Iterate twice during which sample points are refined. 
    % This is done to find the nonsmooth parts of the function
    % and sample more densely there. 
    for it = 1:2 
    num = round(NN/length(xk));             
    Z = [];
    for ii = 1:length(xk)-1
        Z = [Z linspace(xk(ii),xk(ii+1),num)];  % equispaced sampling between
                                                % each pair of reference pts
    end
    Z = unique(Z); Z = Z(:); F = feval(f,Z);
    [r,~,~,~,xk] = aaamn_lawson(F,Z,m,n); % Do AAA-Lawson with updated sample pts
    xk = findReference(f,r,m,n,xk); 
    end    
end    

% Cumulative distribution function-based initialization of the rational
% version of the exchange algorithm. If we want to compute a degree (m,n)
% approximation, the following steps are performed:
%
% 1) compute a type (m-diff,n-diff) best approximation, where 0 < diff < min(m,n);
%
% 2) look at how the references for this approximation are distributed
% inside the domain and determine a piecewise linear approximation of
% their distribution function;
%
% 3) use this function to retrieve a starting reference for a degree
% (m-diff+k,n-diff+k) approximation (k is a divisor of diff, i.e., diff=k*m)
% to f.

% These steps are repeated m times, each time the numerator and denominator
% degrees increasing by k. The step parameter determines the value of k.

function xk = cdfInit(f,m,n,symFlag,opts,step)
    
    stepSize = step;
    if(symFlag > 0) % dealing with symmetry (even or odd)
        stepSize = 2*step;
    end
    text = ['Trying CDF-based initialization with step size ', num2str(stepSize),'...'];
    disp(text);
    
    if(abs(m-n) <= 2)
        % Choose small starting degree
        minValue = min(m,n);
        minValue = minValue - rem(minValue,stepSize);
        k = minValue/stepSize;
        k = k - (3 - step);
        startM = m - stepSize * k;
        startN = n - stepSize * k;
    
        % need an initialization strategy that has a high chance
        % of working without problem for the small degree case;
        % CF is used for now
    
        xk = cfInit(f, startM, startN);
        [~,~,~,~,status] = remezKernel(f, startM, startN, startM+startN, true, xk, opts, 0);
    
        if(status.success == 1)   
            while(startM < m - stepSize && (status.success == 1))
                startM = startM + stepSize;
                startN = startN + stepSize;
            if opts.displayIter
                disp(['CDF (m, n) = ',num2str(startM),', ',num2str(startN)])
            end                
                xk = refGen(f, status.xk, startM + startN + 2, symFlag);
                [~,~,~,~,status] = remezKernel(f, startM, startN, startM+startN, true, xk, opts, 0);
            end
        end
    
        if(status.success == 1)
            xk = refGen(f, status.xk, m + n + 2, symFlag);
        else
            text = ['Initialization failed using CDF with step size ', num2str(stepSize)];
            disp(text);
            [~,~,~,~,status] = remez(f,m+n); xk = status.xk;
        end
    else

        if(m < n)
            minValue = m - rem(m,stepSize);
            k = minValue/stepSize;
            hM = m - stepSize * k;
            hN = n - stepSize * k;
            
            
            hminValue = hN - rem(hN,stepSize);
            hk = hminValue/stepSize;
            startN = hN - stepSize * hk
            xk = cfInit(f, hM, startN);
            [~,~,~,~,status] = remezKernel(f, hM, startN, hM+startN, true, xk, opts, 0);
            
            if(status.success == 1)   
                while(startN < hN - stepSize && (status.success == 1))
                    startN = startN + stepSize;                    
                    xk = refGen(f, status.xk, hM + startN + 2, 0);
                if opts.displayIter
                    disp(['CDF (m, n) = ',num2str(hM),', ',num2str(startN)])
                end                                    
                    [~,~,~,~,status] = remezKernel(f, hM, startN, hM+startN, true, xk, opts, 0);
                end
            end
            
            if(status.success == 1)
                %xk = refGen(f, status.xk, hM + hN + 2, symFlag);
                xk = refGen(f, status.xk, hM + hN + 2, 0);
            else
                [~,~,~,~,status] = remez(f,hM+hN); xk = status.xk;
            end
            
            status.success = 1;  
            while(hM < m - stepSize && (status.success == 1))
                    hM = hM + stepSize;
                    hN = hN + stepSize;
                    xk = refGen(f, status.xk, hM + hN + 2, symFlag);
                if opts.displayIter
                    disp(['CDF (m, n) = ',num2str(hM),', ',num2str(hN)])
                end                                    
                    [~,~,~,~,status] = remezKernel(f, hM, hN, hM+hN, true, xk, opts, 0);
            end
    
            if(status.success == 1)
                %xk = refGen(f, status.xk, m + n + 2, symFlag);
                xk = refGen(f, status.xk, m + n + 2, 0);
            else
                text = ['Initialization failed using CDF with step size ', num2str(stepSize)];
                disp(text);
                [~,~,~,~,status] = remez(f,m+n); xk = status.xk;
            end     
            
        else % m > n
            minValue = n - rem(n,stepSize);
            k = minValue/stepSize;
            startM = m - stepSize * k;
            startN = n - stepSize * k;
            xk = cfInit(f, startM, startN);
            [~,~,~,~,status] = remezKernel(f, startM, startN, startM+startN, true, xk, opts, 0);
    
            if(status.success == 1)   
                while(startM < m - stepSize && (status.success == 1))
                    startM = startM + stepSize;
                    startN = startN + stepSize;
                if opts.displayIter
                    disp(['CDF (m, n) = ',num2str(startM),', ',num2str(startN)])
                end                                    
                    xk = refGen(f, status.xk, startM + startN + 2, symFlag);
                    [~,~,~,~,status] = remezKernel(f, startM, startN, ...
                        startM+startN, true, xk, opts, 0);
                end
            end
    
            if(status.success == 1)
                xk = refGen(f, status.xk, m + n + 2, symFlag);
            else
                text = ['Initialization failed using CDF with step size ', num2str(stepSize)];
                disp(text);
                [~,~,~,~,status] = remez(f,m+n); xk = status.xk;
            end
        end
    end
    
end

function varargout = remezKernel(f, m, n, N, rationalMode, xk, opts, dialogFlag)

% This core function should only ever be called with a nonempty initial set
% of xk reference values
normf = opts.normf;
dom = opts.dom;

% If m=-1, this means f=odd and input (m,n)=(0,n); return constant 0. 
if ( m == -1 )
    q = chebfun(1, dom);
    p = chebfun(0, dom);
    varargout = {p, q, @(x) feval(p, x)./feval(q, x), norm(f,inf), []};    
    return
end

% With zero denominator degree, the denominator polynomial is trivial.
if ( n == 0 )
    q = chebfun(1, dom);
end

% Initial values for some parameters.
iter = 0;                 % Iteration count.
delta = max(normf, eps);  % Value for stopping criterion.
deltamin = inf;           % Minimum error encountered.
diffx = 1;                % Maximum correction to trial reference

xo = xk;

% Print header for text output display if requested.
if ( opts.displayIter && dialogFlag)
    disp('It.   Max(|Error|)     |ErrorRef|    Delta ErrorRef    Delta Ref      m  n')
end

h = -1;
err = normf;
interpSuccess = 1;
% Run the main algorithm.
while ( (abs(abs(h)-abs(err))/abs(err) > opts.tol) && ...
    (iter < opts.maxIter) && (diffx > 0) && (interpSuccess == 1))
    hpre = h;
    if (abs(abs(h)-abs(err))/normf < 1e-14)
        break
    end
    % Compute trial function and levelled reference error.
    if ( n == 0 )
        fk = feval(f, xk);     % Evaluate on the exchange set.
        w = baryWeights(xk);   % Barycentric weights for exchange set.
        [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, dom);
         
        % Perturb exactly-zero values of the levelled error.
        if ( h == 0 )
            h = 1e-19;
        end
 
        rh = @(x) 0;
        % Update the exchange set using the Remez algorithm with full exchange.
        [xk, err, err_handle, ~] = exchange(xk, h, 2, f, p, rh, N + 2, n);
 
        % If overshoot, recompute with one-point exchange.
        if ( err/normf > 1e5 )
            [xk, err, err_handle, ~] = exchange(xo, h, 1, f, p, rh, N + 2, n);
        end
 
        % Update max. correction to trial reference and stopping criterion.
        diffx = max(abs(xo - xk));
        delta = err - abs(h);
 
        % Store approximation with minimum norm.
        if ( delta < deltamin )
            pmin = p;
            errmin = err;
            xkmin = xk;
            deltamin = delta;
        end
         
    else
        err = inf;
        [p, q, rh, h, interpSuccess, ~] = ...
            computeTrialFunctionRational(f, xk, m, n, hpre, dialogFlag);
   
        % Perturb exactly-zero values of the levelled error.
        if ( h == 0 )
            h = 1e-19;
        end
         
        if(interpSuccess == 1)
            [xk, err, err_handle, ~] = exchange(xk, h, 2, f, p, rh, N+2, n);
            diffx = max(abs(xo - xk));
            delta = err - abs(h);
            
            if opts.tol*norm(err_handle(xk),inf) < normf*1e-14
                % relative tolerance is below machine precision, make it
                % reasonable
                opts.tol = normf*1e-13/norm(err_handle(xk),inf);
                opts.tol = min( opts.tol, 0.1 );
            end
        end
    end
 
    % Display diagnostic information as requested.
    if ( opts.plotIter && interpSuccess && dialogFlag)
        doPlotIter(xo, xk, err_handle, h, dom);
    end
 
    if ( opts.displayIter && dialogFlag)
        doDisplayIter(iter, err, h, delta, normf, diffx, m, n);
    end
 
    xo = xk;
    iter = iter + 1;
end
 
if (n == 0)
    % Take best results of all the iterations we ran.
    p = pmin;
    err = errmin;
    xk = xkmin;
    delta = deltamin;
end
 
% Warn the user if we failed to converge.
if ( abs(abs(h)-abs(err))/abs(err) > opts.tol && ...
        abs(abs(h)-abs(err))/normf >= 1e-14 && dialogFlag && interpSuccess)
    warning('CHEBFUN:CHEBFUN:remez:convergence', ...
        ['Remez algorithm did not converge after ', num2str(iter), ...
         ' iterations to the tolerance ', num2str(opts.tol), '.']);
end
 
% Form the outputs.
status.delta = delta/normf;
status.iter = iter;
status.diffx = diffx;
status.xk = xk;
status.success = interpSuccess;
 
if(not(isempty(p)))
    p = simplify(p);
end
if ( rationalMode )
    varargout = {p, q, rh, err, status};
else
    varargout = {p, err, status};
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [m, n, N, rationalMode, symFlag, xk, opts] = parseInputs(f, varargin)

% Detect polynomial / rational approximation type and parse degrees.
if ( ~mod(nargin, 2) ) % Even number of inputs --> polynomial case.
    m = varargin{1};
    n = 0;
    rationalMode = false;
    symFlag = 0;
    varargin = varargin(2:end);
else                   % Odd number of inputs --> rational case.
    [m, n, symFlag] = adjustDegreesForSymmetries(f, varargin{1}, varargin{2});
    if ( n == 0 )
        rationalMode = false;
    else
        rationalMode = true;
    end
    varargin = varargin(3:end);
end

N = m + n;

% Parse name-value option pairs.
if rationalMode
    opts.tol = 1e-4;                        % Relative tolerance for deciding convergence.
    opts.maxIter = 10+round(max(m,n)/2);    % Maximum number of allowable iterations.
else
    opts.tol = 1e-14*(N^2 + 10); % Polynomial case is much more robust. 
    opts.maxIter = 20;           % Maximum number of allowable iterations.
end

opts.displayIter = false;    % Print output after each iteration.
opts.plotIter = false;       % Plot approximation at each iteration.
opts.dom = f.domain([1, end]);
opts.normf = norm(f);
xk = [];

for k = 1:2:length(varargin)
    if ( strcmpi('tol', varargin{k}) )
        opts.tol = varargin{k+1};
    elseif ( strcmpi('maxiter', varargin{k}) )
        opts.maxIter = varargin{k+1};
    elseif ( strcmpi('display', varargin{k}) )
        opts.displayIter = true;
    elseif ( strcmpi('plotfcns', varargin{k}) )
        opts.plotIter = true;
    elseif ( strcmpi('init', varargin{k}) )
        xk = varargin{k+1};
    else
        error('CHEBFUN:CHEBFUN:remez:badInput', ...
            'Unrecognized sequence of input parameters.')
    end
end

end

function [m, n, symFlag] = adjustDegreesForSymmetries(f, m, n)
%ADJUSTDEGREESFORSYMMETRIES   Adjust rational approximation degrees to account
%   for function symmetries.
%
%   [M, N] = ADJUSTDEGREESFORSYMMETRIES(F, M, N) returns new degrees M and N to
%   correct the defect of the rational approximation if the target function is
%   even or odd.  In either case, the Walsh table is covered with blocks of
%   size 2x2, e.g.  for even function the best rational approximant is the same
%   for types [m/n], [m+1/n], [m/n+1] and [m+1/n+1], with m and n even. This
%   strategy is similar to the one proposed by van Deun and Trefethen for CF
%   approximation in Chebfun (see @chebfun/cf.m).

% Sample piecewise-smooth CHEBFUNs.
if ( (numel(f.funs) > 1) || (length(f) > 128) )
  f = chebfun(f, f.domain([1, end]), 128);
end

% Assume no symmetry at the outset.
symFlag = 0;
% Compute the Chebyshev coefficients.
c = chebcoeffs(f, length(f));
c(1) = 2*c(1);

% Check for symmetries and reduce degrees accordingly.
if ( max(abs(c(2:2:end)))/vscale(f) < eps )   % f is even.
    symFlag = 1;
    if ( mod(m, 2) == 1 )
        m = max(m - 1, 0);
    end
    if ( mod(n, 2) == 1 )
        n = max(n - 1, 0);
    end
elseif ( max(abs(c(1:2:end)))/vscale(f) < eps ) % f is odd.
    symFlag = 2;
    if ( ~mod(m, 2) )
        m = m - 1;
    end
    if ( mod(n, 2) )
        n = n - 1;
    end
end

end

function [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, dom)

% Vector of alternating signs.
sigma = ones(N + 2, 1);
sigma(2:2:end) = -1;

h = (w'*fk) / (w'*sigma);                          % Levelled reference error.
pk = (fk - h*sigma);                               % Vals. of r*q in reference.

% Trial polynomial.
p = chebfun(@(x) bary(x, pk, xk, w), dom, m + 1);

end

function [p, q, rh, h, interpSuccess,xsupport] = ...
    computeTrialFunctionRational(f, xk, m, n, hpre, dialogFlag)
% computeTrialFunctionRational finds a rational approximation to f at an 
% iteration of the Remez algorithm. This uses the barycentric representation
% for improved numerical stability. 
% f: function 
% xk: approximate reference points
% m,n: degree

% The function values at the current reference points
fk = feval(f, xk);
% Take barycentric support points to be alternating values of two reference points
xsupport = (xk(1:2:end-1)+xk(2:2:end))/2;  
xadd = (xk(2:2:end-1)+xk(3:2:end))/2;   % when m~=n, we need more support points

if ismember(f.domain(1),xk) == 0        % if endpoints aren't included, add them
    xadd = [(f.domain(1)+xk(1))/2;xadd]; 
end
if ismember(f.domain(end),xk) == 0
    xadd = [(f.domain(end)+xk(end))/2;xadd];
end
num = abs((max(m,n)+1-length(xsupport)));
[xadd, ~] = leja(xadd, 1, num);          % take Leja points from the remaining ref pts

if m~=n
    xsupport = [xsupport;xadd(1:max(m,n)+1-length(xsupport))]; % add any lacking supp pts
end
xsupport = sort(xsupport,'ascend');

if m~=n % force coefficients to lie in null space of Vandermonde
    mndiff = abs(m-n);
    Q = ones(length(xsupport),1);        
    Q = Q/norm(Q);     
    for ii = 2:mndiff
        Qtmp = diag(xsupport)*Q(:,end);
        Qtmp = Qtmp-Q*(Q'*Qtmp);        % orthogonalize
        Qtmp = Qtmp/norm(Qtmp);         % normalize
        Q = [Q Qtmp];
    end
    [Q,~] = qr(Q);
    Qmn = Q(:,mndiff+1:end);            % This is the orthogonal subspace
    [Qmnall,~] = qr(Qmn);
end

    C = zeros(length(xk),length(xsupport)); % Cauchy (basis) matrix
for ii = 1:length(xk)
    C(ii,:) = 1./(xk(ii)-xsupport);
end

% find Delta diag matrix 
wt = zeros( 1,length(xk) );
for ii = 1:length(xk)    
% wt(ii) = prod(xk(ii)-xsupport);
% wxdiff(ii) = prod(xk(ii)-xk([1:ii-1 ii+1:end]));    
% Delta = diag(-(wt.^2)./wxdiff);
% do above in a way that avoids underflow, overflow
    wt(ii) = exp(sum(log(abs(prod(xk(ii)-xsupport)))));
%wxdiff(ii) = exp(sum(log(abs(xk(ii)-xk([1:ii-1 ii+1:end])))));    
    Delta(ii) = -exp(2*sum(log(abs(prod(xk(ii)-xsupport)))) ...
        - sum(log(abs(xk(ii)-xk([1:ii-1 ii+1:end])))));
end
Delta = diag(Delta); 

%DD = diag(1./norms(sqrt(abs(Delta))*C)); % scaling, might help stability
DD = eye(size(C,2));

% prepare QR factorizations; these lead to symmetric eigenproblem
if ( m == n )
    [Q,R] = qr(sqrt(abs(Delta))*C,0);
elseif ( m > n )
    [Q,R] = qr(sqrt(abs(Delta))*C*Qmn,0);    
    [Qall,Rall] = qr(sqrt(abs(Delta))*C*Qmnall,0);
else % m<n
    [Q,R] = qr(sqrt(abs(Delta))*C,0);
    [Qpart,Rpart] = qr(sqrt(abs(Delta))*C*Qmn,0);
end

S = diag((-1).^(0:length(xk)-1));
%Q2 = S*Q; for sanity check svd(Q'*Q2) or svd(Qpart'*Q2) when m<n, should be O(eps)

QSQ = Q'*S*diag(fk)*Q;
QSQ = (QSQ+QSQ')/2; % force symmetry as it's supposed to be

% key operation; this forces (F+hsigma)N=D, where N/D is rational approximant. 
% The eigenvector VR containing the coefficients for 
% D(x)= sum_i VR_{i}/(x-xsupport_{i}). 
[VR,d] = eig(-QSQ); % symmetric eigenproblem
beta = R\VR;        % Denominator coefficients in barycentric form

% obtain alpha (the Numerator coefficients in bary form) from beta
if ( m == n )
    alpha = R\(-Q'*diag(fk)*Q*VR);    
elseif ( m > n )
    alpha = Qmnall*((Rall)\((Qall'*diag(-fk)*Q*VR)));
else % m<n
    alpha = (Rpart)\(Qpart'*diag(-fk)*Q*VR);    
end
vt = [alpha;beta];

% conditioning check, might help
% disp([cond(C) cond(sqrt(abs(Delta))*C) cond(sqrt(abs(Delta))) ...
%     cond(C*DD) cond(sqrt(abs(Delta))*C*DD) m n])

% Among the n+1 eigenvalues, only one can be the solution. The correct
% one needs to have no sign changes in the denominator polynomial
% D(x)*node(x), where node(x) = prod(x-xsupport). 

if ( m <= n ) % values of D at xk
    Dvals = C(:,1:n+1)*(DD*vt(m+1+1:end,:)); 
else
    Dvals = C*(Qmn*vt(m+1+1:end,:)); 
end
node = @(z) prod(z-xsupport); % needed to check sign

nodevec = xk;
for ii = 1:length(xk)
    nodevec(ii) = node(xk(ii));   % values of node polynomial
end
% Find position without sign changes in D*node. Ignore ones with too small
% Dvals. 
pos = find(abs(sum(sign(diag(nodevec)*Dvals))) == m+n+2 & sum(abs(Dvals))>1e-4);  

if isempty(pos)                   % unfortunately no solution with same signs.
    if(dialogFlag)
        disp('Trial interpolant too far from optimal.')
    end
    interpSuccess = 0; 
    
    p = []; q = []; rh = []; h = 1e-19;
    return
elseif ( length(pos) > 1 )        % more than one solution with no sign changes.. try something
    [~,ix] = min(abs(hpre)-diag(abs(d(pos,pos))));
    pos = pos(ix);
end    

h = -d(pos, pos);                 % levelled reference error.

% coefficients for barycentric representations
if ( m <= n )
    wD = (DD*vt(m+2:end,pos));
else
    wD = Qmn*vt(m+2:end,pos);    
end
if ( m >= n )
    wN = DD(1:m+1,1:m+1)*vt(1:m+1,pos);
else
    wN = Qmn*vt(1:m+1,pos);    
end

D = @(x) 0; N = @(x) 0;    % form function handle rh = N/D 
for ii = 1:length(xsupport)
   D = @(x) D(x) + wD(ii)./(x-xsupport(ii));
   N = @(x) N(x) + wN(ii)./(x-xsupport(ii));   
end
D = @(x)-D(x); % flip back sign

rh = @(zz) feval(@rr,zz,xsupport,wN,wD); % rational approximant as function handle
interpSuccess = 1; 

% Form chebfuns of p and q (note: could be numerically unstable but
% provided for convenience)
% Find values of node polynomial at Chebyshev points
x = chebpts(m+n+1);
nodex = zeros(length(x),1); for ii = 1:length(x),    nodex(ii) = node(x(ii)); end 
qvals = nodex.*feval(D,x);  % Values of p and q at Chebyshev points
pvals = nodex.*feval(N,x);
 
p = chebfun(pvals); q = chebfun(qvals);
p = simplify(p); q = simplify(q); % or
% p = chebfun(p,m+1); qp = chebfun(q,n+1); 

end


function r = rr(zz,xsupport,wN,wD)        % evaluate r at zz
zv = zz(:);                               % vectorize zz if necessary
CC = 1./bsxfun(@minus,zv,xsupport.');     % Cauchy matrix 
r = -(CC*(wN))./(CC*wD);                  % barycentric approx as vector
r = reshape(r,size(zz));                  % barycentric approx
end


function [xx, pos] = leja(x, startIndex, nPts) 
% put NPTS from X in a Leja sequence
% starting from x(startIndex)
n = length(x);
p = zeros(n,1);
pos = zeros(nPts, 1);
xx = zeros(nPts, 1);
xx(1) = x(startIndex); 
pos(1) = startIndex;

for j = 2:nPts
    % we want to pick the jth point now:
    for i = 1:n
        %p(i) = prod(abs(x(i) - xx(1:j-1)));
        p(i) = sum(log(abs(x(i) - xx(1:j-1)))); % no overflow
    end  
    [~,pos(j)] = max(p);
    xx(j) = x(pos(j));
end

end


function [xk, norme, err_handle, flag] = exchange(xk, h, method, f, p, rh, Npts, n)
%EXCHANGE   Modify an equioscillation reference using the Remez algorithm.
%   EXCHANGE(XK, H, METHOD, F, P, RH, NPTS, N) performs one step of the Remez algorithm
%   for the best rational approximation of the CHEBFUN F of the target function
%   according to the first method (METHOD = 1), i.e., exchanges only one point,
%   or the second method (METHOD = 2), i.e., exchanges all the reference points.
%   XK is a column vector with the reference, H is the levelled error, P is the
%   numerator, and RH is a function handle, NPTS is the required number of
%   alternation points, and N is the denominator degree.
%
%   [XK, NORME, E_HANDLE, FLAG] = EXCHANGE(...) returns the modified reference
%   XK, the supremum norm of the error NORME (included as an output argument,
%   since it is readily computed in EXCHANGE and is used later in REMEZ), a
%   function handle E_HANDLE for the error, and a FLAG indicating whether there
%   were at least N+2 alternating extrema of the error to form the next
%   reference (FLAG = 1) or not (FLAG = 0).

% Compute extrema of the error.
if(n == 0) % polynomial case
    e_num = diff(f) - diff(p);
    rts = roots(e_num, 'nobreaks');
    rr = [f.domain(1) ; rts; f.domain(end)];
    % Function handle output for evaluating the error.
    err_handle = @(x) feval(f, x) - feval(p, x);
else       % rational case
    rr = findExtrema(f, rh, xk);
    err_handle = @(x) feval(f, x) - rh(x);
end

% Select exchange method.
if ( method == 1 )                           % One-point exchange.
    [~, pos] = max(abs(feval(err_handle, rr)));
    pos = pos(1);
else                                           % Full exchange.
    pos = find(abs(err_handle(rr)) >= abs(h)); % Values above levelled error
end

% Add extrema nearest to those which are candidates for exchange to the
% existing exchange set.
[r, m] = sort([rr(pos) ; xk]);
v = ones(Npts, 1);
v(2:2:end) = -1;
er = [feval(err_handle, rr(pos)) ; v*h];
er = er(m);

% Delete repeated points.
repeated = diff(r) == 0;
r(repeated) = [];
er(repeated) = [];

% Determine points and values to be kept for the reference set.
s = r(1);    % Points to be kept.
es = er(1);  % Values to be kept.
for i = 2:length(r)
    if ( (sign(er(i)) == sign(es(end))) && (abs(er(i)) > abs(es(end))) )
        % Given adjacent points with the same sign, keep one with largest value.
        s(end) = r(i);
        es(end) = er(i);
    elseif ( sign(er(i)) ~= sign(es(end)) )
        % Keep points which alternate in sign.
        s = [s ; r(i)];    %#ok<AGROW>
        es = [es ; er(i)]; %#ok<AGROW>
    end
end



% Of the points we kept, choose n + 2 consecutive ones that include the maximum
% of the error.
[norme, index] = max(abs(es));
d = max(index - Npts + 1, 1);
if ( Npts <= length(s) )
    xk = s(d:d+Npts-1);
    flag = 1;
else
    xk = s;
    flag = 0;
end

end

function rts = findExtrema(f,rh,xk)
% finds all the local maxima and minima
% of f-p./q
% xk is the present reference
% rh is a handle to p/q
% h is the current leveled error

err_handle = @(x) feval(f, x) - rh(x);
sample_points = linspace(f.domain(1),f.domain(end),1000);
scale_of_error = norm(err_handle(sample_points),inf);
relTol =  1e-15 * (vscale(f)/scale_of_error);  % adjust to make this right!
rts = [];

doms = unique([f.domain(1); xk; f.domain(end)]).';

% Initial trial
if ( isempty(xk) )
    ek = chebfun(@(x) err_handle(x), f.domain, 'eps', relTol, 'splitting', 'on');
    rts = roots(diff(ek), 'nobreaks');
else
    for k = 1:length(doms)-1
        ek = chebfun(@(x) err_handle(x), [doms(k), doms(k+1)], 33, 'eps', 1e-12); 
        ek = simplify(ek);
        rts = [rts; roots(diff(ek), 'nobreaks')];  %#ok<AGROW>
    end    
end

% Append end points of the domain:
rts = unique([f.domain(1) ; rts; f.domain(end)]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find extrema of error function in AAA-Lawson
function xk = findReference(f,r,m,n,z) 
    % f: function  
    % r: rational approximant 
    % m,n: (m,n) is type of rational approximant
    % z: barycentric support points 
    % OUTPUT: reference points xk
        
    xk = findExtrema(f,r, sort(z,'ascend'));      % find extrema points as usual
    
    % Deal with length(xk) not equal to the desired m+n+2
    if length(xk) > m+n+2 % reduce reference pts becuse too many found 
   
        xkdiff = diff(xk);                        
        [~,ix] = sort(xkdiff,'descend'); % take those with largest gaps
        xk = [xk(1);xk(1+ix(1:m+n+1))];
        xk = sort(xk,'ascend');        
    elseif length(xk) < m+n+2 % increase reference pts becuse too few found 

        xkdiff = diff(xk);                        
        add = m+n+2-length(xk);          % we need to add this many reference points 
        [~,ix] = sort(xkdiff,'descend'); % take those with largest gaps and fill midpoints
        xk = [xk;(xk(ix(1:add))+xk(ix(1:add)+1))/2];
        xk = sort(xk,'ascend');
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for displaying diagnostic information.

% Function called when opts.plotIter is set.
function doPlotIter(xo, xk, err_handle, h, dom)

xxk = linspace(dom(1), dom(end), 10000);
plot(xo, err_handle(xo), 'or', 'MarkerSize', 4)   % Old reference.
holdState = ishold;
hold on
plot(xk, err_handle(xk), '*k', 'MarkerSize', 4)   % New reference.
plot(xxk, err_handle(xxk))               % Error function.
plot(xxk, ones(size(xxk))*h,'r');
plot(xxk, -ones(size(xxk))*h,'r');
if ( ~holdState )                        % Return to previous hold state.
    hold off
end
xlim(dom)
legend('Current Ref.', 'Next Ref.', 'Error')
drawnow

end

% Function called when opts.displayIter is set.
function doDisplayIter(iter, err, h, delta, normf, diffx, m, n)

disp([num2str(iter,'%3g'), '      ', num2str(err, '%5.4e'), '      ', ...
    num2str(abs(h), '%5.4e'), '      ', ...
    num2str(delta/normf, '%5.4e'), '      ', num2str(diffx, '%5.4e'),...
    '    ', num2str(m, '%4g'),'  ', num2str(n, '%4g')])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions used by the CDF-based initialization routine

% Constructs a piecewise linear function which interpolates the data
% (xd(i), yd(i)) when i goes from 1 to length(xd). It then computes
% the value of this piecewise linear approximation at the xi nodes
% (ni is the number of xi nodes)
function yi = pwiselin(xd, yd, ni, xi)

  nd = length(xd);
  xd = xd(:);
  yd = yd(:);
  xi = xi(:);

  if ( nd == 1 )
    yi(1:ni,1) = yd;
    return
  end

  [~, ~,k] = histcounts(xi, xd);

  k ( k == 0 ) = 1;
  k ( k == nd ) = nd - 1;

  t = ( xi - xd(k,1) ) ./ ( xd(k+1,1) - xd(k,1) );
  yi = ( 1 - t ) .* yd(k) + t .* yd(k+1);
  
  return
end

% Generate a set of n reference points which follow a distribution
% of the xk nodes. The symType flag tells us if we are dealing with
% an even or odd function, in which case the new reference nodes are
% taken to be symmetric with respect to the middle of the approximation
% domain
function nxk = refGen(f, xk, n, symType)

xx = linspace(-1,1,length(xk));

if(symType == 0)

    nxk = pwiselin(xx, xk, n, linspace(-1,1,n));

% handling of even symmetries

elseif(symType == 1)

    halfSize = length(xx)/2;
    halfn = n/2;


    if (xk(1) == f.domain(1))
        nxk = pwiselin(xx(halfSize+1:end), xk(halfSize+1:end),halfn, ...
            linspace(xx(halfSize+1),xx(end),halfn));
        nxk = [nxk; -nxk(2:end); f.domain(1)];
        nxk = sort(nxk,'ascend');
    elseif (xk(end) == f.domain(end))
        nxk = pwiselin(xx(1:halfSize), xk(1:halfSize),halfn, ...
            linspace(xx(1),xx(halfSize),halfn));
        nxk = [nxk; -nxk(1:end-1); f.domain(end)];
        nxk = sort(nxk,'ascend');
    else
        nxk = pwiselin(xx, xk, n, linspace(-1,1,n));
    end

% handling of odd symmetries
else

halfSize = (length(xx)-1) / 2;
halfn = (n-1)/2;

    if (xk(1) == f.domain(1))
        nxk = pwiselin(xx(halfSize+2:end), xk(halfSize+2:end),halfn, ...
            linspace(xx(halfSize+2),xx(end),halfn));
        nxk = [nxk; -nxk(1:end); f.domain(1)];
        nxk = sort(nxk,'ascend');
    elseif (xk(end) == f.domain(end))
        nxk = pwiselin(xx(1:halfSize+1), xk(1:halfSize+1),halfn, ...
            linspace(xx(1),xx(halfSize+1),halfn));
        nxk = [nxk; -nxk(1:end); f.domain(end)];
        nxk = sort(nxk,'ascend');
    else
        nxk = pwiselin(xx, xk, n, linspace(-1,1,n));
    end

end

end
