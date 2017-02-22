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
%   [...] = REMEZ(..., 'init', XK) allows the user to specify the vector XK as
%   the starting reference.
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
%   [1] B. Beckermann, S. Filip and Y. Nakatsukasa, manuscript in preparation.
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
[m, n, N, rationalMode, polyOutput, symFlag, xk, opts] = parseInputs(f, varargin{:});

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
        if(polyOutput)
            varargout = {p,err,status};
        else
            varargout = {p,q,p,err,status};
        end
    else % rational case
        % A first attempt is using CF as an initial guess
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
        q = chebfun('1');
        if(polyOutput)
            varargout = {p,err,status};
        else
            varargout = {p,q,p,err,status};
        end
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
    warning off
    if ( numel(f.funs) == 1 )
        [p, q] = cf(f, m, n);
    else
        [p, q] = cf(f, m, n, 50*(m+n) );
    end
    warning on
    pqh = @(x) feval(p, x)./feval(q, x);
    [xk, ~, ~, flag] = exchange([], 0, 2, f, p, pqh, m+n+2, n);

    % If the above procedure failed to produce a reference
    % with enough oscillation points, use polynomial Remez.
    if ( flag == 0 )
        [~,~,status] = remez(f,m+n); xk = status.xk;
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
% version of the exchange algorithm.
function xk = cdfInit(f,m,n,symFlag,opts,step)
    newlineCounter = 0;
    stepSize = step; % Increase in degree of numerator and/or denominator
                     % at each step.
    if(symFlag > 0) % Dealing with symmetry (even or odd).
        stepSize = 2*step;
    end
    text = ['Trying CDF-based initialization with step size ', ...
                                num2str(stepSize),'...'];
    disp(text);
    
    % The approximation is close to being diagonal; start from an
    % approximation where both the numerator and denominator degrees
    % are decreased by the same value.
    if(abs(m-n) <= 2)

        minValue = min(m,n);
        minValue = minValue - rem(minValue,stepSize);
        k = minValue/stepSize;
        k = k - (3 - step);
        % Starting small degree problem.
        startM = m - stepSize * k;
        startN = n - stepSize * k;
    
        % We need an initialization strategy that has a high chance
        % of working without problem for the small degree case;
        % CF is used for now.
        xk = cfInit(f, startM, startN);
        [~,~,~,~,status] = remezKernel(f, startM, startN, startM+startN, ...
                                        true, xk, opts, 0);
        % If Remez worked on the small degree problem, start increasing
        % the degrees in both the numerator and denominator.
        if(status.success == 1)   
            while(startM < m - stepSize && (status.success == 1))
                startM = startM + stepSize;
                startN = startN + stepSize;
                newlineCounter = newlineCounter + 1;
                if(newlineCounter == 10)
                    newlineCounter = 0;
                    fprintf('\n');
                end
                fprintf('(%d,%d) ',startM, startN);
                % Use the distribution information from the previous
                % reference to construct a starting reference for the new,
                % larger degree problem.
                xk = refGen(f, status.xk, startM + startN + 2, symFlag);
                [~,~,~,~,status] = remezKernel(f, startM, startN, ...
                                        startM+startN, true, xk, opts, 0);
            end
        end
    

        if(status.success == 1)
            xk = refGen(f, status.xk, m + n + 2, symFlag);
            fprintf('(%d,%d)\n',m, n);
        else
            % There was a failure somewhere, use reference from polynomial
            % Remez (might work sometimes).
            fprintf('\n');
            text = ['Initialization failed using CDF with step size ', ...
                                            num2str(stepSize)];
            disp(text);
            [~,~,status] = remez(f,m+n); xk = status.xk;
        end
    else
        % Similar strategy to the diagonal case.
        if(m < n)
            % Construct the 'corner' instance.
            % (m - stepSize*k, n - stepSize*k), where m - stepSize * k will
            % usually be 0 or 1.
            minValue = m - rem(m,stepSize);
            k = minValue/stepSize;
            hM = m - stepSize * k;
            hN = n - stepSize * k;
            
            % Now decrease the degree only in the denominator.
            hminValue = hN - rem(hN,stepSize);
            hk = hminValue/stepSize;
            startN = hN - stepSize * hk;
            xk = cfInit(f, hM, startN);
            [~,~,~,~,status] = remezKernel(f, hM, startN, hM+startN, ...
                                            true, xk, opts, 0);
            % Construct approximations by successively increasing the
            % denominator degree.
            if(status.success == 1)   
                while(startN < hN - stepSize && (status.success == 1))
                    startN = startN + stepSize;                    
                    xk = refGen(f, status.xk, hM + startN + 2, 0);
                    newlineCounter = newlineCounter + 1;
                    if(newlineCounter == 10)
                        newlineCounter = 0;
                        fprintf('\n');
                    end
                    fprintf('(%d,%d) ',hM, startN);
                    [~,~,~,~,status] = remezKernel(f, hM, startN, ...
                                            hM+startN, true, xk, opts, 0);
                end
            end
            
            if(status.success == 1)
                %xk = refGen(f, status.xk, hM + hN + 2, symFlag);
                xk = refGen(f, status.xk, hM + hN + 2, 0);
            else
                [~,~,status] = remez(f,hM+hN); xk = status.xk;
            end
            
            % Go to the initial degree by now simultaneously increasing
            % both numerator and denominator degree.
            status.success = 1;  
            while(hM < m - stepSize && (status.success == 1))
                    hM = hM + stepSize;
                    hN = hN + stepSize;
                    xk = refGen(f, status.xk, hM + hN + 2, symFlag);
                    newlineCounter = newlineCounter + 1;
                    if(newlineCounter == 10)
                        newlineCounter = 0;
                        fprintf('\n');
                    end
                    fprintf('(%d,%d) ',hM, hN);
                    [~,~,~,~,status] = remezKernel(f, hM, hN, hM+hN, ...
                                                true, xk, opts, 0);
            end
    
            if(status.success == 1)
                %xk = refGen(f, status.xk, m + n + 2, symFlag);
                xk = refGen(f, status.xk, m + n + 2, 0);
                fprintf('(%d,%d)\n',m, n);
            else
                fprintf('\n');
                text = ['Initialization failed using CDF with step size ', ...
                                num2str(stepSize)];
                disp(text);
                [~,~,status] = remez(f,m+n); xk = status.xk;
            end     
            
        else % m > n
            % Construction of the 'corner' instance by decreasing the
            % degree in both numerator and denominator.
            minValue = n - rem(n,stepSize);
            k = minValue/stepSize;
            startM = m - stepSize * k;
            startN = n - stepSize * k;
            xk = cfInit(f, startM, startN);
            [~,~,~,~,status] = remezKernel(f, startM, startN, ...
                                       startM+startN, true, xk, opts, 0);
    
            if(status.success == 1)   
                while(startM < m - stepSize && (status.success == 1))
                    startM = startM + stepSize;
                    startN = startN + stepSize;
                    newlineCounter = newlineCounter + 1;
                    if(newlineCounter == 10)
                        newlineCounter = 0;
                        fprintf('\n');
                    end
                    fprintf('(%d,%d) ',startM, startN);
                    xk = refGen(f, status.xk, startM + startN + 2, symFlag);
                    [~,~,~,~,status] = remezKernel(f, startM, startN, ...
                                        startM+startN, true, xk, opts, 0);
                end
            end
    
            if(status.success == 1)
                xk = refGen(f, status.xk, m + n + 2, symFlag);
                fprintf('(%d,%d)\n',m, n);
            else
                fprintf('\n');
                text = ['Initialization failed using CDF with step size ', ...
                          num2str(stepSize)];
                disp(text);
                [~,~,status] = remez(f,m+n); xk = status.xk;
            end
        end
    end
    fprintf('\n');
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
    % Approximation error is at the level of machine precision, stop.
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
                % Relative tolerance below machine precision, make it
                % reasonable.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [m, n, N, rationalMode, polyOutput, symFlag, xk, opts] = parseInputs(f, varargin)

% Detect polynomial / rational approximation type and parse degrees.
polyOutput = true;
if ( ~mod(nargin, 2) ) % Even number of inputs --> polynomial case.
    m = varargin{1};
    n = 0;
    rationalMode = false;
    symFlag = 0;
    varargin = varargin(2:end);
else                   % Odd number of inputs --> rational case.
    polyOutput = false;
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
    opts.tol = 1e-4;                     % Relative tolerance for deciding convergence.
    opts.maxIter = 10+round(max(m,n)/2); % Maximum number of allowable iterations.
else
    opts.tol = 1e-14*(N^2 + 10); % Polynomial case is much more robust. 
    opts.maxIter = 30;           % Maximum number of allowable iterations.
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

h = (w'*fk) / (w'*sigma);           % Levelled reference error.
pk = (fk - h*sigma);                % Vals. of r*q in reference.

% Trial polynomial.
p = chebfun(@(x) bary(x, pk, xk, w), dom, m + 1);

end

function [p, q, rh, h, interpSuccess,xsupport] = ...
    computeTrialFunctionRational(f, xk, m, n, hpre, dialogFlag)
% computeTrialFunctionRational finds a rational approximation to f at an 
% iteration of the Remez algorithm. It uses a barycentric representation
% for improved numerical stability. 
% f: function 
% xk: approximate reference points
% m,n: degree

% The function values at the current reference points
fk = feval(f, xk);
% Take barycentric support points to be alternating values of two reference points
xsupport = (xk(1:2:end-1)+xk(2:2:end))/2;  
xadd = (xk(2:2:end-1)+xk(3:2:end))/2; % when m~=n, we need more support points

if ismember(f.domain(1),xk) == 0      % if endpoints aren't included, add them
    xadd = [(f.domain(1)+xk(1))/2;xadd]; 
end
if ismember(f.domain(end),xk) == 0
    xadd = [(f.domain(end)+xk(end))/2;xadd];
end
num = abs((max(m,n)+1-length(xsupport)));
[xadd, ~] = leja(xadd, 1, num);  % take Leja points from the remaining ref pts

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

if isempty(pos)  % Unfortunately, no solution with same signs.
    if(dialogFlag)
        disp('Trial interpolant too far from optimal.')
    end
    interpSuccess = 0; 
    
    p = []; q = []; rh = []; h = 1e-19;
    return
elseif ( length(pos) > 1 ) % more than one solution with no sign changes...
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

% Form chebfuns of p and q (note: could be numerically unstable, but
% provided for convenience)
% Find values of node polynomial at Chebyshev points
x = chebpts(m+n+1,f.domain([1,end]));
nodex = zeros(length(x),1); for ii = 1:length(x),    nodex(ii) = node(x(ii)); end 
qvals = nodex.*feval(D,x);  % Values of p and q at Chebyshev points
pvals = nodex.*feval(N,x);
 
p = chebfun(pvals,f.domain([1,end])); q = chebfun(qvals,f.domain([1,end]));
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
else       % Rational case.
    rr = findExtrema(f, rh, xk);
    err_handle = @(x) feval(f, x) - rh(x);
end

% Select exchange method.
if ( method == 1 )                           % One-point exchange.
    [~, pos] = max(abs(feval(err_handle, rr)));
    pos = pos(1);
else                                           % Full exchange.
    pos = find(abs(err_handle(rr)) >= abs(h)); % Values above levelled error.
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



% Of the points we kept, choose n + 2 consecutive ones that include the
% maximum of the error.
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
% Finds all the local maxima and minima of f-rh.
% xk is the present reference
% rh is a handle to p/q

err_handle = @(x) feval(f, x) - rh(x);
sample_points = linspace(f.domain(1),f.domain(end),5000);
scale_of_error = norm(err_handle(sample_points),inf);
relTol =  1e-15 * (vscale(f)/scale_of_error);
rts = [];

doms = unique([f.domain(1); xk; f.domain(end)]).';

% Initial trial
if ( isempty(xk) )
    warning off
    ek = chebfun(@(x) err_handle(x), f.domain, 'eps', relTol, 'splitting', 'on');
    warning on
    rts = roots(diff(ek), 'nobreaks');
else
    for k = 1:length(doms)-1
        ek = chebfun(@(x) err_handle(x), [doms(k), doms(k+1)], 33, 'eps', 1e-12); 
        ek = simplify(ek);
        rts = [rts; roots(diff(ek), 'nobreaks')];  %#ok<AGROW>
    end    
end

% Append end points of the domain.
rts = unique([f.domain(1) ; rts; f.domain(end)]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find extrema of error function in AAA-Lawson
function xk = findReference(f,r,m,n,z) 
    % f: function  
    % r: rational approximant 
    % m,n: (m,n) is the type of our rational approximant
    % z: barycentric support points 
    % OUTPUT: reference points xk
        
    xk = findExtrema(f,r, sort(z,'ascend')); % Find extrema points as usual.
    
    % Deal with length(xk) not equal to the desired m+n+2.
    if length(xk) > m+n+2 % Reduce reference pts becuse too many found. 
   
        xkdiff = diff(xk);                        
        [~,ix] = sort(xkdiff,'descend'); % Take those with largest gaps.
        xk = [xk(1);xk(1+ix(1:m+n+1))];
        xk = sort(xk,'ascend');        
    elseif length(xk) < m+n+2 % Increase # of reference pts if too few found. 

        xkdiff = diff(xk);                        
        add = m+n+2-length(xk);  % We need to add this many reference points. 
        % Take those with largest gaps and fill midpoints.
        [~,ix] = sort(xkdiff,'descend');
        xk = [xk;(xk(ix(1:add))+xk(ix(1:add)+1))/2];
        xk = sort(xk,'ascend');
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



function [r, pol, res, zer, z, Z, f, w, wf, errvec, p, q] = aaamn_lawson(F, varargin)
%AAAMN_Lawson   near-best rational approximation of F. 
% 
% R = aaamn_lawson(F) computes a rational aproximant of (default) type
% (10,10) on the default interval [-1,1]
%
% [R, POL, RES, ZER] = aaamn_lawson(...) outputs the poles, residues and zeros
% of R. The poles, zeros will approximate those of F (not well if R-F is not small)
% 
% [R, POL, RES, ZER, z, Z, F, W,WF, ERRVEC, P, Q] = aaamn_lawson(...) 
% outputs additionally the sample points Z, support points z, values
% F=f(Z), weights w and wf (see below) and AAA errvec, and p,q = chebfuns 
% s.t. r= p/q (note this can be numerically unstable)
%
% [...] = aaamn_lawson(F,m,n) specifies the type to (m,n).
%
% [...] = aaamn_lawson(F,Z,m,n) also specifies the sample points Z
% (recommended). 
%
% [...] = aaamn_lawson(F,Z,m,n,'plot','on') will plot the error functions
% as the Lawson iterations proceed. 
%
% [...] = aaamn_lawson(F,m,n,'dom',[-1,2]) specifies the domain (this has no effect 
% if Z is specified)
%
% [...] = aaamn_lawson(F,m,n,'tol',1e-5) specifies the Lawson iterate
% relative stopping criterion (here to 1e-5)
%
% [...] = aaamn_lawson(F,m,n,'iter',10) limits the Lawson maximum
% iterations to 10. 
% 
% 
% The algorithm first finds a AAA rational approximant to F \approx r(Z), then
% attempts to refine the approximant by a Lawson process, i.e. an iterative
% reweighting. 
% 
% This code is designed for computing good reference points for the
% rational remez code to follow, but can be used independently for
% constructing a rational approximation r that can be much closer than AAA
% to the best rational approximant. 
%
% Input:  Z = vector of sample points
%         F = vector of data values, or a function handle
%         m, n: max type is (m,n), set to 10 if omitted
%         tol = relative tolerance tol, default: 1e-13 
%         Lawsoniter: max. iteration number of Lawson updates (default 10)
%         doplot: 1 to plot error curve history (default 0)
%     R = AAAMN_Lawson(F, Z, m, n, NAME, VALUE) sets the following parameters:
%   - 'tol', TOL: relative tolerance for Lawson iteration (default 1e-5)
%   - 'iter', IT: maximal number of Lawson iterations (default MMAX = max([5 min([20, m, n])])).
%         
%
% Output: r = AAA-Lawson approximant to F (function handle)
%         pol,res,zer = vectors of poles, residues, zeros
%         errvec = vector of errors at each step
%         z,f,w,wf = vectors of support pts, function values, weights
%         s.t. r = N/D, N(x) = sum_i wf(i)/(x-z(i)) and 
%         D(x) = sum_i w(i)/(x-z(i)).
%         p,q = chebfuns s.t. r= p/q (note this can be numerically unstable)
%
% Examples:
%    r = aaamn_lawson(@abs,10,10)
%    r = aaamn_lawson(@abs,10,10,'plot','on')
%    [r, pol, res, zer, z, f, w, wf, errvec, p, q] = aaamn_lawson(@abs,10,10,'plot','on','dom',[-1 2])
%    r = aaamn_lawson(@exp,4,2,'plot','on','dom',[-1 2])
%    r = aaamn_lawson(@(x)log(1.1-x),5,5,'plot','on')
%
%    f = chebfun(@(x)-1./(log(abs(x)).^2),[-.1,.1],'splitting','on'); 
%    [r,pol,res] = aaamn_lawson(f,linspace(-.1,.1,1e4),18,18,'plot','on')
%
%

% parse inputs
[F, Z, m, n, Lawsoniter, tolLawson, doplot, tol ] = ...
    parseInputslawson(F, varargin{:});

M = length(Z);                            % number of sample points
mmax = m+1; nmax = n+1;                   % for coding convenience
if ( (nargin < 6) || isempty(Lawsoniter) )% number of Lawson updates
    Lawsoniter = max([5 min([20,mmax,nmax])]); 
end 
if ~isfloat(F), F = feval(F,Z); end        % convert function handle to vector
 Z = Z(:); F = F(:);                       % work with column vectors
 SF = spdiags(F,0,M,M);                    % left scaling matrix
 J = 1:M;                                  % indices that are not support pts
 z = []; f = []; C = [];                   % initializations
 errvec = []; R = mean(F); 
for mn = 1:max(mmax,nmax)
  [~,j] = max(abs(F-R));                  % select next support point
  z = [z; Z(j)];                          % update set of support points
  f = [f; F(j)];                          % update set of data values
  J(J==j) = [];                           % update index vector
  C = [C 1./(Z-Z(j))];                    % next column of Cauchy matrix
  Sf = diag(f);                           % right scaling matrix
  A = SF*C - C*Sf;                        % Loewner matrix
    
        if ( mn > min(nmax,mmax) ) % nondiagonal case, find projection subspace 
             if mmax < nmax
             q = f(:);
             else
             q = ones(length(z),1);
             end
             Q = orthspace(z,mn-min(mmax,nmax),q);   % projection subspace 
             [~,~,V] = svd(A(J,:)*Q,0);              % SVD on projected subspace
             w = Q*V(:,end);
        else             
             [~,~,V] = svd(A(J,:),0);               % SVD, no projection needed
             w = V(:,mn);                           % weight vector             
        end     
  wf = w.*f;
  N = C*(w.*f); D = C*w;                  % numerator and denominator
  R = F; R(J) = N(J)./D(J);               % rational approximation
  err = norm(F-R,inf);
  errvec = [errvec; err];                 % max error at sample points
  if ( err < tol*norm(F,inf) ), break, end    % stop if converged
end
    r = @(zz) feval(@rrint,zz,z,w,f);            % AAA approximant as function handle
    Rori = R;
    
% now start Lawson, in this mode we leave interpolation and work with
% 'alpha-beta' mode. 
          wei = ones(length(J),1);
          nrmbest = inf;
          if ( mn > min(nmax,mmax) )  % Deal with projection for m neq n                                             
              if mn>nmax
                A =[SF*C*Q -C];        
              else % need to redefine Q as not the same as AAA above        
                q = ones(length(z),1);
                Q = orthspace(z,mn-min(mmax,nmax),q);     % projection subspace                              
                A =[SF*C -C*Q];                   
              end
          else
            A =[SF*C -C];   % diagonal case
          end
          
          rate = 1;         % default Lawson rate, will shrink if not converging
          nrmincreased = 0; % initialization    
	      for it = 1:Lawsoniter
              weiold = wei; 
              wei = wei .* power(abs(F(J)-R(J)),rate); % update Lawson weights
              wei = wei/sum(wei);                      % normalize 
              if norm(weiold-wei)/norm(wei)< tolLawson % declare Lawson converged
                  break
              end
              D = spdiags(sqrt(wei),0,length(wei),length(wei)); % diagonal weight matrix

              [~,~,V] = svd(D*A(J,:),0);     % weighted least-squares via SVD
              
          if ( mn > min(nmax,mmax) )    % deal with nondiagonal case
              if ( mn > nmax )
                w = Q*V(1:nmax,end); wf = V(nmax+1:end,end);            
              else
                w = V(1:nmax,end); wf = Q*V(nmax+1:end,end); 
              end
          else
            w = V(1:mn,end); wf = V(mn+1:2*mn,end);            
          end                      
            f = wf./w;                          % for compatibility with interpolatory-aaa                        
            N = C*wf; D = C*w;                  % numerator and denominator               
            R = F; R(J) = N(J)./D(J);           % rational approximation
            err = norm(F-R,inf);            
            errvec = [errvec; err];             % max error at sample points                            
            if ( err < nrmbest )    % adopt best so far
            nrmbest = norm(F-R,'inf'); 
            weibest = wei;                      % store best weight
            r = @(zz) feval(@rrab,zz,z,w,wf,f); % AAA approximant as function handle
            else
                nrmincreased = nrmincreased + 1;
            end
            if ( nrmincreased >= 3 )     % perhaps not converging,
            rate = max( rate/2,0.01 );   % make Lawson update conservative
            if doplot
                warning(['Lawson rate made conservative to ',num2str(rate)])
            end
              nrmincreased = 0; 
            end

            if doplot  % plot error functions (hopefully near-equioscillating)
                subplot(2,1,1)
                plot(Z,F-Rori,'r.','markersize',8)
                title('AAA error')
                grid on, hold on
                if exist('hh','var')
                    set(hh,'color',(0.8)*[1 1 1]); 
                end
                subplot(2,1,2)
                title('AAA-Lawson error')
                hh = plot(Z,F-R,'k.','markersize',8);                        
                grid on, hold on
                h2 = plot(z,0*z,'m.','markersize',12);
                ylim(err*[-1 1]); drawnow, shg            
                if ( it == Lawsoniter ) % plot best function 
                plot(Z,F-r(Z),'b.','markersize',10);
                end
            end            
	      end          

    % compute poles and roots
    B = eye(mn+1); B(1,1) = 0;                 
    E = [0 wf.'; ones(mn,1) diag(z)];      
    zer = eig(E,B); zer = zer(~isinf(zer));   % zeros  
    E = [0 w.'; ones(mn,1) diag(z)];      
    pol = eig(E,B); pol = pol(~isinf(pol));   % poles
    dz = 1e-5*exp(2i*pi*(1:4)/4);
    res = r(bsxfun(@plus,pol,dz))*dz.'/4;     % residues
    
    if ( nargout > 8 )  % form p and q via N/D, NOTE: not always stable
    D = @(x) 0;    N = @(x) 0;     % form r=N/D in barycentric form
    for ii = 1:length(z)
       D = @(x) D(x) + w(ii)./(x-z(ii));
       N = @(x) N(x) + wf(ii)./(x-z(ii));               
    end
    dom = [min(real(Z)) max(real(Z))];  % set domain
    x = chebpts(mmax+nmax+1,dom);
    node = @(x) prod(x-z);              % needed to check sign
    nodex = zeros(length(x),1);         % setup node values    
    for ii = 1:length(x), nodex(ii) = node(x(ii)); end
    qvals = nodex.*feval(D,x);          % values of p,q
    pvals = nodex.*feval(N,x);
    p = chebfun(pvals,dom); q = chebfun(qvals,dom); % form chebfuns    
    end
end    


% parse Inputs for Lawson:

function [F, Z, m, n, Lawsoniter, tolLawson, doplot, tol ] = ...
    parseInputslawson(F, varargin)
% Input parsing for AAAmn_lawson.

% Check if F is empty:
if ( isempty(F) )
    error('CHEBFUN:aaamn_lawson:emptyF', 'No function given.')
elseif ( isa(F, 'chebfun') )
    if ( size(F, 2) ~= 1 )
        error('CHEBFUN:aaamn_lawson:nColF', 'Input chebfun must have one column.')
    end
end

% Domain:
if ( isa(F, 'chebfun') )
    dom = F.domain([1, end]);
else
    dom = [-1, 1];
end

% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    if length(varargin{1})>2   % sample points Z given. 
    Z = varargin{1};
    varargin(1) = [];    
    else                       % sample points not given.
        
    end
end

% m,n
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    if length(varargin{1})>1 % input (f,Z,[m n])
    mn = varargin{1};
    m = mn(1); n = mn(2);
    varargin(1) = [];    
    elseif isfloat(varargin{2}) % input (f,Z,m,n)
    m = varargin{1};        
    n = varargin{2};    
    varargin([1, 2]) = [];
    else
    m = varargin{1};        
    n = m; % input (f,Z,m), default to diagonal type (m,m)
    varargin(1) = [];
    end
end

if ( ~exist('m', 'var') ) 
     warning('CHEBFUN:aaamn_lawson: type (m,n) not specified, default to (10,10)')
     m = 10; n = 10; 
end

% Set defaults for other parameters:
tolLawson = 1e-5;                       % Relative tolerance for Lawson update.
tol = 1e-15;                            % AAA tolerance
Lawsoniter = max([5 min([20, m, n])]);  % Maximum number of terms.
doplot = 0;                             % Don't plot intermediate functions unless specified

% Check if parameters have been provided:
while ( ~isempty(varargin) )
    if ( strncmpi(varargin{1}, 'tol', 3) |  strncmpi(varargin{1}, 'tolLawson', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tolLawson = varargin{2};   % Lawson tolerance
        end
        varargin([1, 2]) = [];

    elseif ( strncmpi(varargin{1}, 'iter', 4) | strncmpi(varargin{1}, 'maxit', 5))
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            Lawsoniter = varargin{2};  % maximum Lawson iterations
        else
        warning(['CHEBFUN:aaamn_lawson:iter unspecified, use default itermax ', num2str(Lawsoniter)])
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'dom', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 2]) )
            dom = varargin{2};
        end
        varargin([1, 2]) = [];
        if ( isa(F, 'chebfun') )
            if ( ~isequal(dom, F.domain([1, end])) )
                warning('CHEBFUN:aaamn_lawson:dom', ...
                    ['Given domain does not match the domain of the chebfun.\n', ...
                    'Results may be inaccurate.'])
            end
        end
        
    elseif strncmpi(varargin{1}, 'plot', 4)  % plot error functions
        if isfloat(varargin{2})
            doplot = varargin{2};
        elseif ( strncmpi(varargin{2}, 'true', 4) | strncmpi(varargin{2}, 'on', 2) )
            doplot = 1;
        end
        varargin([1, 2]) = [];                
    else
        error('CHEBFUN:aaamn_lawson:UnknownArg', 'Argument unknown.')
    end
end

% Deal with Z and F:
if ( ~exist('Z', 'var') && isfloat(F) )
    % F is given as data values, pick same number of sample points:
    Z = linspace(dom(1), dom(2), length(F)).';
end

if ( exist('Z', 'var') )
    % Work with column vector:
    Z = Z(:);
    M = length(Z);
    
    % Function values:
    if ( isa(F, 'function_handle') || isa(F, 'chebfun') )
        % Sample F on Z:
        F = F(Z);
    elseif ( isnumeric(F) )
        % Work with column vector and check that it has correct length.
        F = F(:);
        if ( length(F) ~= M )
            error('CHEBFUN:aaamn_lawson:lengthFZ', ...
                'Inputs F and Z must have the same length.')
        end
    else
        error('CHEBFUN:aaamn_lawson:UnknownF', 'Input for F not recognized.')
    end
    
else
    % Z was not given.  Set flag that Z needs to be determined.
    % Also set Z and M since they are needed as output.
    % in AAA this is done adaptively. This can be done with Lawson, but
    % probably safe to take as many points as reasonably possible here. 
    Z = linspace(dom(1), dom(end), 4000).';    
end

end % End of PARSEINPUT().


%% generate function handle, interpolatory mode
function r = rrint(zz,z,w,f)                 % evaluate r at zz
zv = zz(:);                               % vectorize zz if necessary
CC = 1./bsxfun(@minus,zv,z.');            % Cauchy matrix 
r = (CC*(w.*f))./(CC*w);                  % AAA approx as vector
ii = find(isnan(r));                      % Find values NaN = 0/0 if any
for j = 1:length(ii)
  r(ii(j)) = f(find(zv(ii(j))==z));       % Force interpolation there
end
r = reshape(r,size(zz));                  % AAA approx
end

% generate function handle, non-interpolatory mode
function r = rrab(zz,z,w,wf,f)                 % evaluate r at zz
zv = zz(:);                               % vectorize zz if necessary
CC = 1./bsxfun(@minus,zv,z.');            % Cauchy matrix 
r = (CC*(wf))./(CC*w);                  % AAA approx as vector

ii = find(isnan(r));                      % Find values NaN = 0/0 if any
for j = 1:length(ii)
  r(ii(j)) = f(find(zv(ii(j))==z));       % Force interpolation there
end

r = reshape(r,size(zz));                  % AAA approx
end

% find null space for nondiagonal case m~=n
function Q = orthspace(z,dim,q)       % orthonormal projection space for (m,n)
if ( dim == 0 ), Q = eye(length(z)); end 
if ( nargin < 3 ), q = ones(length(z),1); end
                Q = q/norm(q);
                for ii = 2:dim               % orthogonal complement via Lanczos-type process
                Qtmp = diag(z)*Q(:,end);
                Qtmp = Qtmp - Q*(Q'*Qtmp);   % orthogonalize
                Qtmp = Qtmp/norm(Qtmp);      % normalize
                Q = [Q Qtmp]; 
                end
            [Q,~] = qr(Q); Q = conj(Q(:,dim+1:end)); % desired null space
end

