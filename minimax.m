function varargout = minimax(f, varargin)
%MINIMAX   Best polynomial or rational approximation for real valued
%          continuous functions. This code supersedes REMEZ. 
%
%   P = MINIMAX(F, M) computes the minimax polynomial approximation of
%   degree M to the real function F using the Remez algorithm. F can
%   be either a CHEBFUN, a function handle or a string representation
%   of the function to approximate.  P is a CHEBFUN.
%
%   [P, Q] = MINIMAX(F, M, N) computes the minimax rational approximation
%   P/Q of type (M, N).  P and Q are CHEBFUNs, but in difficult cases
%   working with P/Q is numerically unstable.
% 
%   [P, Q, R_HANDLE] = MINIMAX(F, M, N) additionally returns a numerically
%   stable function handle for evaluating P/Q (based on a barycentric
%   representation).
%
%   [...] = MINIMAX(..., [A, B]) takes the approximation domain to be
%   [A, B]. If a domain is not specified and F is a CHEBFUN, then the
%   domain of F is used. In all other cases, [-1, 1] is used.
%
%   [...] = MINIMAX(..., 'tol', TOL) uses the value TOL as the termination
%   tolerance on the relative equioscillation error.  The default is 
%   approximately 1e-12 for polynomial approximation and 1e-4 for
%   rational approximation.
%
%   [...] = MINIMAX(..., 'display', 'iter') displays output at each
%   iteration.
%
%   [...] = MINIMAX(..., 'maxiter', MAXITER) sets the maximum number of
%   allowable iterations to MAXITER.
%
%   [...] = MINIMAX(..., 'init', XK) allows the user to specify the vector
%   XK as the starting reference.
%
%   [...] = MINIMAX(..., 'plot', 'on'), or equivalently 
%   [...] = MINIMAX(..., 'plotfcns', 'error') plots the error after each
%   iteration while the algorithm executes.
%
%   [...] = MINIMAX(..., 'silent') turns off all messages regarding the
%   execution of the algorithm.
%
%   [P, ERR] = MINIMAX(...) and [P, Q, R_HANDLE, ERR] = MINIMAX(...)
%   return the maximum error estimate ERR.
%
%   [P, ERR, STATUS] = MINIMAX(...) and [P, Q, R_HANDLE, ERR, STATUS] =
%   MINIMAX(...) return a structure array STATUS with the following fields:
%       STATUS.DELTA - Tolerance obtained.
%       STATUS.ITER  - Number of iterations performed.
%       STATUS.DIFFX - Maximum correction in last trial reference.
%       STATUS.XK    - Last trial reference on which the error
%                      equioscillates.
%   In case we are doing rational approximation (denominator degree >=1),
%   two extra fields are computed:
%       STATUS.POL   - Poles of the minimax approximation.
%       STATUS.ZER   - Zeros of the minimax approximation.
%
%   This code is highly reliable for polynomial approximation but may
%   sometimes have difficulties in the rational case, though we believe
%   it is the most powerful rational minimax code available.
%
% Examples:
%   x = chebfun('x'); f = abs(x);
%   p = minimax(f, 20); plot(f-p)
%   [p, q, rh] = minimax(f, 10, 10);
%   xx = linspace(-1,1,10000); plot(xx, f(xx)-rh(xx))
%
% References:
%
%   [1] B. Beckermann, S. Filip, Y. Nakatsukasa and L. N. Trefethen,
%   "Rational minimax approximation via adaptive barycentric
%   representations", arXiv:1705.10132.
%
%   [2] R. Pachon and L. N. Trefethen, "Barycentric-Remez algorithms for
%   best polynomial approximation in the chebfun system", BIT Numerical
%   Mathematics, 49:721-742, 2009.
%
% See also AAA, CF, CHEBPADE, PADEAPPROX, RATINTERP, POLYFIT, POLYFITL1.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(f,'chebfun') ) % check if input is chebfun; if not, look for
                         % splitting points
    dom = [];
    domIndex = 0;
    
    for k = 1:length(varargin) % look for domain
        if isfloat(varargin{k})
            if length(varargin{k}) == 2
                dom = varargin{k};
                domIndex = k;
            end
        end
    end
    
    if isempty(dom) % domain not provided, default to [-1,1]
        dom = [-1 1];
    else
        varargin(domIndex) = [];
    end
    if ( ischar(f) )
        fHandle = str2op(vectorize(f));
    else
        fHandle = f;
    end
    f = chebfun(f, dom, 'splitting', 'on');
else % f is a chebfun input
    fHandle = @(x) feval(f, x);
end

if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:minimax:real', ...
        'MINIMAX only supports real valued functions.');
end

if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:minimax:quasi', ...
        'MINIMAX does not currently support quasimatrices.');
end

% Parse the inputs.
[m, n, N, rationalMode, polyOutput, symFlag, xk, opts] = ...
                                   parseInputs(f, varargin{:});

% If m = -1, this means f = odd and input (m,n) = (0,n); return constant 0. 
if ( m == -1 )
    q = chebfun(1, f.domain([1, end]));
    p = chebfun(0, f.domain([1, end]));
    varargout = {p, q, @(x) feval(p, x)./feval(q, x), norm(f,'inf'), []};    
    return
end

if ( isempty(xk) ) % no initial reference is given by the user
    % Several initialization attempts are made
    if ( n == 0 ) % polynomial case
        % Try Chebyshev points
        xk = chebpts(N + 2, f.domain([1, end]));
        [p,err,status] = minimaxKernel(f, fHandle,m, n, N, rationalMode,...
                                       xk, opts, 1);
        q = chebfun('1', f.domain([1, end]));
        if ( polyOutput )
            varargout = {p, err, status};
        else
            varargout = {p, q, p, err, status};
        end
    else % rational case
        % A first attempt is using CF as an initial guess                
        try
            xk = cfInit(f, fHandle, m, n);
            [p,q,rh,err,status] = minimaxKernel(f, fHandle,m, n, N, ...
                                                rationalMode, xk, opts, 1);
            varargout = {p, q, rh, err, status};
        catch
            if ~opts.silentFlag
                disp(['CF-based initialization failed,' ...
                    ' turning to AAA-Lawson...']);
            end
            status.success = 0; % CF didn't work
        end        
        
        % If CF doesn't give a satisfactory answer, we try AAA-Lawson
        if ~status.success
            if ~opts.silentFlag
                disp('Trying AAA-Lawson-based initialization...');
            end
            xk = AAALawsonInit(f, fHandle, m, n);
            [p,q,rh,err,status] = minimaxKernel(f, fHandle, m, n, N, ...
                                                rationalMode, xk, opts, 1);
            varargout = {p, q, rh, err, status};
        end
        

        % A final attempt using cumulative distribution functions
        if ~status.success
            xk = cdfInit(f, fHandle, m, n, symFlag, opts, 1);
            [p,q,rh,err,status] = minimaxKernel(f, fHandle, m, n, N, ...
                                                rationalMode, xk, opts, 1);
            varargout = {p, q, rh, err, status};            
        end
        

        
        if ~status.success
            xk = cdfInit(f, fHandle, m, n, symFlag, opts, 2);
            [p,q,rh,err,status] = minimaxKernel(f, fHandle, m, n, N, ...
                                                rationalMode, xk, opts, 1);
            varargout = {p, q, rh, err, status};
        end
        
        if ~status.success % all attempts failed
        	error('CHEBFUN:CHEBFUN:minimax:failure', ...
               ['MINIMAX failed to produce the best approximant. ' ...
                'If the accuracy is close to machine precision, ' ...
                'it may be that what you''ve asked for is unachievable ' ...
                'in floating-point arithmetic. In such a case, try ' ...
                'reducing the degree to get a clean best approximant.'])
        end
    end
else  % the user has also given a starting reference
    if ( n == 0 )
        [p,err,status] = minimaxKernel(f, fHandle, m, n, N, ...
                                       rationalMode, xk, opts, 1);
        q = chebfun('1', f.domain([1,end]));
        if ( polyOutput )
            varargout = {p, err, status};
        else
            varargout = {p, q, p, err, status};
        end
    else
        [p,q,rh,err,status] = minimaxKernel(f, fHandle, m, n, N, ...
                                            rationalMode, xk, opts, 1);
        varargout = {p, q, rh, err, status};
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implementing the different initialization strategies

% CF-based initialization
function xk = cfInit(f, fHandle, m, n)
    warning off
    if ( numel(f.funs) == 1 )
        [p, q] = cf(f, m, n);
        pqh = @(x) feval(p, x)./feval(q, x);
        [xk, ~, ~, flag] = exchange([], 0, 2, f, fHandle, p, pqh, m+n+2, n);
    else
        try
            [p, q] = cf(f, m, n, 50*(m+n) );
            pqh = @(x) feval(p, x)./feval(q, x);
            [xk, ~, ~, flag] = exchange([], 0, 2, f, fHandle, p, pqh, ...
                                        m+n+2, n);
        catch ME  % an error occured when calling cf (ignore it)
            flag = 0;
        end
    end
    warning on

    % If the above procedure failed to produce a reference
    % with enough oscillation points, use polynomial Remez.
    if ( flag == 0 )
        [~,~,status] = minimax(f, m+n); xk = status.xk;
    end
end

% Now turn to initialization via AAA-Lawson, this is more expensive than CF
% but less so than CDF (which follows if this fails). 
% 
function xk = AAALawsonInit(f,fHandle, m,n) % AAA-Lawson initialization for
                                            % functions with breakpoints
    NN = max(10*max(m,n),round(1e5/max(m,n))); 
    dom = domain(f);
    Z = linspace(dom(1), dom(end), NN); 
    F = fHandle(Z);
    [r,~,~,~,xk] = aaamn_lawson(F, Z, m, n);    % 1st AAA-Lawson    
    xk = findReference(f, fHandle, r, m, n, xk);
    
    % Iterate twice during which sample points are refined. 
    % This is done to find the nonsmooth parts of the function
    % and sample more densely there. 
    for it = 1:2 % (maybe helps to run more)
    num = round(NN/length(xk));             
    Z = [];
    for ii = 1:length(xk)-1
        Z = [Z linspace(xk(ii), xk(ii+1),num)]; % equispaced sampling
                                                % between each pair of
                                                % reference pts
    end
    Z = unique(Z); Z = Z(:); F = feval(f, Z);
    [r,~,~,~,xk] = aaamn_lawson(F, Z, m, n); % Do AAA-Lawson with updated
                                             % sample pts
    xk = findReference(f, fHandle, r, m, n, xk); 
    end    
end    

% Cumulative distribution function-based initialization of the rational
% version of the exchange algorithm.
function xk = cdfInit(f, fHandle, m, n, symFlag, opts, step)
    newlineCounter = 0;
    stepSize = step; % Increase in degree of numerator and/or denominator
                     % at each step.
    if ( symFlag > 0 ) % Dealing with symmetry (even or odd).
        stepSize = 2*step;
    end
    if ~opts.silentFlag
        text = ['Trying CDF-based initialization with step size ', ...
                                                num2str(stepSize),'...'];
        disp(text);
    end
    
    % The approximation is close to being diagonal; start from an
    % approximation where both the numerator and denominator degrees
    % are decreased by the same value.
    if ( abs(m-n) <= 2 )
        minValue = min(m, n);
        minValue = minValue - rem(minValue, stepSize);
        k = minValue / stepSize;
        k = k - (3 - step);
        % Starting small degree problem.
        startM = m - stepSize * k;
        startN = n - stepSize * k;
    
        % We need an initialization strategy that has a high chance
        % of working without problem for the small degree case;
        % CF is used for now.
        xk = cfInit(f, fHandle, startM, startN);
        [~,~,~,~,status] = minimaxKernel(f, fHandle, startM, startN, ...
                                        startM+startN, true, xk, opts, 0);
        % If Remez worked on the small degree problem, start increasing
        % the degrees in both the numerator and denominator.
        if status.success   
            while(startM < m - stepSize && (status.success == 1))
                startM = startM + stepSize;
                startN = startN + stepSize;
                newlineCounter = newlineCounter + 1;
                if(newlineCounter == 10)
                    newlineCounter = 0;
                    if ~opts.silentFlag
                        fprintf('\n');
                    end
                end
                if ~opts.silentFlag
                    fprintf('(%d,%d) ',startM, startN);
                end
                % Use the distribution information from the previous
                % reference to construct a starting reference for the new,
                % larger degree problem.
                xk = refGen(f, status.xk, startM + startN + 2, symFlag);
                [~,~,~,~,status] = minimaxKernel(f, fHandle, startM, ...
                                    startN, startM+startN, true, xk, ...
                                    opts, 0);
            end
        end
    

        if status.success
            xk = refGen(f, status.xk, m + n + 2, symFlag);
            if ~opts.silentFlag
                fprintf('(%d,%d)\n',m, n);
            end
        else
            % There was a failure somewhere, use reference from polynomial
            % Remez (might work sometimes).
            if ~opts.silentFlag
                fprintf('\n');
                text = ['Initialization failed using CDF with step', ...
                                    ' size ', num2str(stepSize) '...'];
                disp(text);
            end
            [~,~,status] = minimax(f, m+n); xk = status.xk;
        end
    else
        % Similar strategy to the diagonal case.
        if ( m < n )
            % Construct the 'corner' instance.
            % (m - stepSize*k, n - stepSize*k), where m - stepSize * k will
            % usually be 0 or 1.
            minValue = m - rem(m, stepSize);
            k = minValue / stepSize;
            hM = m - stepSize * k;
            hN = n - stepSize * k;
            
            % Now decrease the degree only in the denominator.
            hminValue = hN - rem(hN, stepSize);
            hk = hminValue / stepSize;
            startN = hN - stepSize * hk;
            xk = cfInit(f, fHandle, hM, startN);
            [~,~,~,~,status] = minimaxKernel(f, fHandle, hM, startN, ...
                                        hM+startN, true, xk, opts, 0);
            % Construct approximations by successively increasing the
            % denominator degree.
            if status.success   
                while ( startN < hN - stepSize && status.success )
                    startN = startN + stepSize;                    
                    xk = refGen(f, status.xk, hM + startN + 2, 0);
                    newlineCounter = newlineCounter + 1;
                    if ( newlineCounter == 10 )
                        newlineCounter = 0;
                        if ~opts.silentFlag
                            fprintf('\n');
                        end
                    end
                    if ~opts.silentFlag
                        fprintf('(%d,%d) ',hM, startN);
                    end
                    [~,~,~,~,status] = minimaxKernel(f, fHandle, hM, ...
                                        startN, hM+startN, true, xk, ...
                                        opts, 0);
                end
            end
            
            if ~status.success
                [~,~,status] = minimax(f, hM + hN);
            end
            
            % Go to the initial degree by now simultaneously increasing
            % both numerator and denominator degree.
            status.success = 1;  
            while ( hM < m - stepSize && status.success )
                    hM = hM + stepSize;
                    hN = hN + stepSize;
                    xk = refGen(f, status.xk, hM + hN + 2, symFlag);
                    newlineCounter = newlineCounter + 1;
                    if(newlineCounter == 10)
                        newlineCounter = 0;
                        if ~opts.silentFlag
                            fprintf('\n');
                        end
                    end
                    if ~opts.silentFlag
                        fprintf('(%d,%d) ', hM, hN);
                    end
                    [~,~,~,~,status] = minimaxKernel(f, fHandle, hM, ...
                                                    hN, hM+hN, true, ...
                                                    xk, opts, 0);
            end
    
            if status.success
                xk = refGen(f, status.xk, m + n + 2, 0);
                if ~opts.silentFlag
                    fprintf('(%d,%d)\n', m, n);
                end
            else
                if ~opts.silentFlag
                    fprintf('\n');
                    text = ['Initialization failed using CDF with ', ...
                                'step size ', num2str(stepSize)];
                    disp(text);
                end
                [~,~,status] = minimax(f, m+n); xk = status.xk;
            end     
            
        else % m > n
            % Construction of the 'corner' instance by decreasing the
            % degree in both numerator and denominator.
            minValue = n - rem(n, stepSize);
            k = minValue / stepSize;
            startM = m - stepSize * k;
            startN = n - stepSize * k;
            xk = cfInit(f, fHandle, startM, startN);
            [~,~,~,~,status] = minimaxKernel(f, fHandle, startM, ...
                                    startN, startM + startN, true, ...
                                    xk, opts, 0);
    
            if status.success   
                while ( startM < m - stepSize && status.success )
                    startM = startM + stepSize;
                    startN = startN + stepSize;
                    newlineCounter = newlineCounter + 1;
                    if ( newlineCounter == 10 )
                        newlineCounter = 0;
                        if ~opts.silentFlag
                            fprintf('\n');
                        end
                    end
                    if ~opts.silentFlag
                        fprintf('(%d,%d) ', startM, startN);
                    end
                    xk = refGen(f, status.xk, startM + startN + 2, ...
                                    symFlag);
                    [~,~,~,~,status] = minimaxKernel(f, fHandle, ...
                                        startM, startN, startM+startN, ...
                                        true, xk, opts, 0);
                end
            end
    
            if status.success
                xk = refGen(f, status.xk, m + n + 2, symFlag);
                if ~opts.silentFlag
                    fprintf('(%d,%d)\n', m, n);
                end
            else
                if ~opts.silentFlag
                    fprintf('\n');
                    text = ['Initialization failed using CDF with ', ...
                            'step size ', num2str(stepSize)];
                    disp(text);
                end
                [~,~,status] = minimax(f, m+n); xk = status.xk;
            end
        end
    end
    if ~opts.silentFlag
        fprintf('\n');
    end
end

function varargout = minimaxKernel(f, fHandle, m, n, N, rationalMode, ...
                                   xk, opts, dialogFlag)

% This core function should only ever be called with a nonempty initial set
% of xk reference values
normf = opts.normf;
dom = opts.dom;

% If m = -1, this means f = odd and input (m,n) = (0,n); return constant 0. 
if ( m == -1 && dialogFlag)
    q = chebfun(1, dom);
    p = chebfun(0, dom);
    varargout = {p, q, @(x) feval(p, x)./feval(q, x), norm(f,inf), []};    
    return
elseif ( m == -1 )
    q = [];
    p = [];
    varargout = {p, q, [], [], []};    
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
    disp('It.   Max(|Error|)     |ErrorRef|    Delta ErrorRef    Delta Ref     m  n')
end

err = normf;
% Initialise the levelled error such that one iteration always executes
h = 2*err + 1;
interpSuccess = 1;
% Run the main algorithm.
while ( (abs(abs(h)-abs(err))/abs(err) > opts.tol) && ...
    (iter < opts.maxIter) && (diffx > 0) && interpSuccess )
    hpre = h;
    % Approximation error is at the level of machine precision, stop.
    if ( abs(abs(h)-abs(err))/normf < 1e-14 )
        break
    end
    % Compute trial function and levelled reference error.
    if ( n == 0 )
        fk = fHandle(xk);      % Evaluate on the exchange set.
        w = baryWeights(xk);   % Barycentric weights for exchange set.
        [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, dom);
         
        % Perturb exactly-zero values of the levelled error.
        if ( h == 0 )
            h = 1e-19;
        end
 
        rh = @(x) 0;
        % Update the exchange set using the Remez algorithm
        % with full exchange rule.
        [xk, err, err_handle, ~] = exchange(xk, h, 2, f, fHandle, p, ...
                                            rh, N + 2, n);
 
        % If overshoot, recompute with one-point exchange rule.
        if ( err/normf > 1e5 )
            [xk, err, err_handle, ~] = exchange(xo, h, 1, f, fHandle, ...
                                            p, rh, N + 2, n);
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
        [p, q, rh, h, interpSuccess, tk, alpha, beta] = ...
            computeTrialFunctionRational(f, fHandle, xk, m, n, hpre, ...
                                         dialogFlag, opts.silentFlag);
   
        % Perturb exactly-zero values of the levelled error.
        if ( h == 0 )
            h = 1e-19;
        end
         
        if(interpSuccess == 1)
            [xk, err, err_handle, ~] = exchange(xk, h, 2, f, fHandle, ...
                                            p, rh, N+2, n);
            diffx = max(abs(xo - xk));
            delta = err - abs(h);
            
            if opts.tol*norm(err_handle(xk),inf) < normf*1e-14
                % Relative tolerance below machine precision, make it
                % reasonable.
                opts.tol = normf*1e-13/norm(err_handle(xk),inf);
                opts.tol = min(opts.tol, 0.1);
            end
        end
    end
 
    % Display diagnostic information as requested.
    if ( opts.plotIter && interpSuccess && dialogFlag )
        doPlotIter(xo, xk, err_handle, h, dom);
    end
 
    if ( opts.displayIter && dialogFlag )
        doDisplayIter(iter, err, h, delta, normf, diffx, m, n);
    end
 
    xo = xk;
    iter = iter + 1;
end
 
if ( n == 0 )
    % Take best results of all the iterations we ran.
    p = pmin;
    err = errmin;
    xk = xkmin;
    delta = deltamin;
end
 
% Warn the user if we failed to converge.
if ( abs(abs(h)-abs(err))/abs(err) > opts.tol && ...
     abs(abs(h)-abs(err))/normf >= 1e-14 && dialogFlag && interpSuccess )
    warning('CHEBFUN:CHEBFUN:minimax:convergence', ...
        ['minimax algorithm did not converge after ', num2str(iter), ...
         ' iterations to the tolerance ', num2str(opts.tol), '.']);
end
 
% Form the outputs.
status.delta = delta/normf;
status.iter = iter;
status.diffx = diffx;
status.xk = xk;
status.success = interpSuccess;
% Compute the poles and zeros in case of a rational approximation
if status.success && dialogFlag && rationalMode
    [status.zer, status.pol] = pzeros(tk, alpha, beta, rh, m, n, dom);
else
    status.zer = []; status.pol = [];
end
 
if( ~isempty(p))
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

function [m, n, N, rationalMode, polyOutput, symFlag, xk, opts] = ...
                                             parseInputs(f, varargin)

opts.silentFlag = 0;
isSilent = 0;
for k = 1:length(varargin)
    if ( ischar(varargin{k}) && strcmpi('silent', varargin{k}) )
        opts.silentFlag = k;
        isSilent = 1;
    end
end

if opts.silentFlag
    varargin(opts.silentFlag) = [];
end

% Detect polynomial / rational approximation type and parse degrees.
polyOutput = true;
if ( ~mod(nargin - isSilent, 2) ) % Even number of inputs --> polynomial
                                  % case.
    m = varargin{1};
    n = 0;
    rationalMode = false;
    symFlag = 0;
    varargin = varargin(2:end);
else                              % Odd number of inputs --> rational case.
    polyOutput = false;
    [m, n, symFlag] = adjustDegreesForSymmetries(f, varargin{1}, ...
                                                 varargin{2});
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
    opts.tol = 1e-4;                     % Relative tolerance for deciding
                                         % convergence.
    opts.maxIter = 10+round(max(m,n)/2); % Maximum number of allowable
                                         % iterations.
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
    elseif ( strcmpi('plot', varargin{k}) )
        opts.plotIter = true;
    elseif ( strcmpi('init', varargin{k}) )
        xk = varargin{k+1};
    else
        error('CHEBFUN:CHEBFUN:minimax:badInput', ...
            'Unrecognized sequence of input parameters.')
    end
end

end

function [m, n, symFlag] = adjustDegreesForSymmetries(f, m, n)
%ADJUSTDEGREESFORSYMMETRIES   Adjust rational approximation degrees to
%   account for function symmetries.
%
%   [M, N] = ADJUSTDEGREESFORSYMMETRIES(F, M, N) returns new degrees M and
%   N to correct the defect of the rational approximation if the target
%   function is even or odd.  In either case, the Walsh table is covered
%   with blocks of size 2x2, e.g.  for even function the best rational
%   approximant is the same for types [m/n], [m+1/n], [m/n+1] and
%   [m+1/n+1], with m and n even. This strategy is similar to the one
%   proposed by van Deun and Trefethen for CF approximation in Chebfun
%   (see @chebfun/cf.m).

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
if ( max(abs(c(2:2:end)))/vscale(f) < eps )     % f is even.
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

function [p, q, rh, h, interpSuccess,xsupport, wN, wD] = ...
    computeTrialFunctionRational(f, fHandle, xk, m, n, hpre, ...
                                 dialogFlag, silentFlag)
% computeTrialFunctionRational finds a rational approximation to f at an 
% iteration of the Remez algorithm. It uses a barycentric representation
% for improved numerical stability. 
% f:          function 
% xk:         approximate reference points
% m,n:        type
% hpre:       levelled error at the previous iteration

% The function values at the current reference points
fk = fHandle(xk);
% Take barycentric support points to be alternating values of two
% reference points
    xsupport = xk(2:2:end);
    xsuppind = 2:2:length(xk);
    xadd = xk(1:2:end);
    xother = xadd; 
    xotherind = 1:2:length(xk);
    
if m~=n % need to add more support points
    xadd = xother;
    % take Leja points from the remaining ref pts and add
    [xadd, ~] = leja(xadd, 1, length(xadd));  
    xsupport = [xsupport;xadd(1:max(m,n)+1-length(xsupport))];   
    xother = zeros(m+n+2-max(m,n)-1,1); 
    xotherind = zeros(m+n+2-max(m,n)-1,1); 
    xsuppind = zeros(max(m,n)+1,1);
    
    iother = 1; isupp = 1;
    for ii = 1:length(xk) % keep indices for later use
        if ~ismember(xk(ii),xsupport)
            xother(iother) = xk(ii);
            xotherind(iother) = ii;
            iother = iother+1;
        else
            xsuppind(isupp) = ii;
            isupp = isupp+1;
        end
    end    
    
end
    xsupport = sort(xsupport,'ascend');
    
if m~=n
    % projection matrices that force coefficients to lie in null space 
    % of Vandermonde matrix
    % projection subspace
    Qmn = orthspace(xsupport,abs(m-n),ones(length(xsupport),1));     
    [Qmnall,~] = qr(Qmn);            
end

    C = 1./bsxfun(@minus,xk,xsupport.');    % Cauchy matrix

    % form matrix Cstar = sqrt(|Delta|)*C        
    % Cstar(ii,jj) = |wt(xi)/sqrt(wx'(xi))|/(xi-tj)    
    Xkdiff = abs(bsxfun(@minus, xk, xk.'));
    Xkdiff(eye(size(Xkdiff))~=0) = 1;       % inf to 0
    Xtdiff = abs(bsxfun(@minus, xother, xsupport.'));
    ST = sum(log(Xtdiff.')); SX = sum(log(Xkdiff));
    VV = exp(ST.'-0.5*SX(xotherind).');
    VV = VV*ones(1,length(xsupport));
    Div = bsxfun(@minus,xother,xsupport.');
    C1 = VV./Div; % odd columns of Cstar

    % diag elements Cstar(ii,jj) = |wt(xi)/sqrt(wx'(xi))|/(xi-tj)    
    Xtdiff = abs(bsxfun(@minus, xsupport, xsupport.'));
    Xtdiff(eye(size(Xtdiff))~=0) = 1;
    ST = sum(log(Xtdiff.'));SX = sum(log(Xkdiff));
    C2 = diag(exp(ST.'-0.5*SX(xsuppind).'));
    Cstar = [C1;C2];
    Cstar(xsuppind,:) = C2;
    Cstar(xotherind,:) = C1;    
    
% prepare QR factorizations; these lead to a symmetric eigenproblem
% we need to be careful how to do QR as rows have
% large dynamical range (though orthogonal columns when m=n)    
% do Householder QR with row sorting, better than [Q,R] = qr(Cstar,0);
if ( m == n )
    nrm = zeros(1, size(Cstar, 1));
    for ii = 1:length(nrm)
        nrm(ii) = norm(Cstar(ii,:));
    end
    [~,ix] = sort(nrm, 'descend');
    %[~,ix] = sort(norms(Cstar'), 'descend');
    [Q,R] = qr(Cstar(ix,:),0);
    ixx(ix) = 1:length(ix);    Q = Q(ixx,:);        
    
    %{
    nrm = zeros(1,m+1); % Cholesky QR with col-scaling, this works too
    for ii = 1:m+1, nrm(ii) = norm(Cstar(:,ii));    end
    Cstar = Cstar/diag(nrm);    % this normalized Cstar is orthogonal
    CTC = Cstar'*Cstar; 
    R = chol(CTC);    Q = Cstar/R;    R = R*diag(nrm);
    %}
elseif ( m > n )
    nrm = zeros(1, size(Cstar, 1));
    for ii = 1:length(nrm)
        nrm(ii) = norm(Cstar(ii,:));
    end
    [~,ix] = sort(nrm, 'descend');
    %[~,ix] = sort(norms(Cstar'),'descend');
    [Q,R] = qr(Cstar(ix,:)*Qmn,0);    
    [Qall,Rall] = qr(Cstar(ix,:)*Qmnall,0);    
    ixx(ix) = 1:length(ix);   
    Q = Q(ixx,:);        Qall = Qall(ixx,:);        
    
else % m<n
    nrm = zeros(1, size(Cstar, 1));
    for ii = 1:length(nrm)
        nrm(ii) = norm(Cstar(ii,:));
    end
    [~,ix] = sort(nrm, 'descend');
    %[~,ix] = sort(norms(Cstar'),'descend');
    [Q,R] = qr(Cstar(ix,:),0);    
    [Qpart,Rpart] = qr(Cstar(ix,:)*Qmn,0);    
    ixx(ix) = 1:length(ix);    
    Q = Q(ixx,:);        Qpart = Qpart(ixx,:);        
end

S = diag((-1).^(0:length(xk)-1));
%Q2 = S*Q; for sanity check svd(Q'*Q2) or svd(Qpart'*Q2) when m<n,
% should be O(eps)

QSQ = Q'*S*diag(fk)*Q;
QSQ = (QSQ+QSQ')/2; % force symmetry as it's supposed to be

% key operation; this forces (F+hsigma)N=D, where N/D is rational
% approximant. The eigenvector VR containing the coefficients for 
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

% Among the n+1 eigenvalues, only one can be the solution. The correct
% one needs to have no sign changes in the denominator polynomial
% D(x)*node(x), where node(x) = prod(x-xsupport). 

if ( m <= n ) % values of D at xk
    % Dvals = C(:,1:n+1)*vt(m+1+1:end,:); 
    bet = vt(m+1+1:end,:);
else
    bet = Qmn*vt(m+1+1:end,:);
end
    Dvals = C*bet; 
    
node = @(z) prod(z-xsupport); % needed to check sign
nodevec = xother;
for ii = 1:length(xother)
    nodevec(ii) = node(xother(ii));   % values of node polynomial
end
% Find position without sign changes in D*node. 
% Evaluate this separately for xsupport and xother.
% Ignore ones with too small Dvals. 

checksign = zeros(length(xk),n+1);
% sign at supp pts
checksign(1:length(xsupport),:) = ...
    diag((-1).^(max(m,n):max(m,n)+length(xsupport)-1))*bet; 
% sign at other pts
signs = sign(diag(nodevec)*Dvals(xotherind,:));
checksign(length(xsupport)+1:end,:) = signs;
pos = find(abs(sum(sign(checksign))) == m+n+2 & sum(abs(Dvals))>1e-7);

if isempty(pos)  % Unfortunately, no solution with same signs.
                 % Try old remez.

% Take barycentric support points to be alternating values of two
% reference points
xsupport = (xk(1:2:end-1)+xk(2:2:end))/2;  
xadd = (xk(2:2:end-1)+xk(3:2:end))/2; % when m~=n, we need more support
                                      % points

if ismember(f.domain(1),xk) == 0      % if endpoints aren't included,
                                      % add them
    xadd = [(f.domain(1)+xk(1))/2;xadd]; 
end
if ismember(f.domain(end),xk) == 0
    xadd = [(f.domain(end)+xk(end))/2;xadd];
end
num = abs((max(m,n)+1-length(xsupport)));
[xadd, ~] = leja(xadd, 1, num);  % take Leja points from the
                                 % remaining ref pts

if m~=n
    % add any lacking supp pts
    xsupport = [xsupport;xadd(1:max(m,n)+1-length(xsupport))];
end
xsupport = sort(xsupport, 'ascend');

C = 1./bsxfun(@minus,xk,xsupport.');

% find Delta diag matrix 
Delta = zeros( 1,length(xk) );
for ii = 1:length(xk)    
% wt(ii) = prod(xk(ii)-xsupport);
% wxdiff(ii) = prod(xk(ii)-xk([1:ii-1 ii+1:end]));    
% Delta = diag(-(wt.^2)./wxdiff); do in a way that avoids
%                                 underflow, overflow
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
% Q2 = S*Q; for sanity check svd(Q'*Q2) or svd(Qpart'*Q2) when m<n,
% should be O(eps)

QSQ = Q'*S*diag(fk)*Q;
QSQ = (QSQ+QSQ')/2; % force symmetry as it's supposed to be

% key operation; this forces (F+hsigma)N=D, where N/D is rational
% approximant. The eigenvector VR for which beta=R\VR contains the
% coefficients for D(x)= sum_i beta_{i}/(x-xsupport_{i}). 
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

if ( m <= n ) % values of D at xk
    Dvals = C(:,1:n+1)*(DD*vt(m+1+1:end,:)); 
else
    Dvals = C*(Qmn*vt(m+1+1:end,:)); 
end
node = @(z) prod(z-xsupport);     % needed to check sign

nodevec = xk;
for ii = 1:length(xk)
    nodevec(ii) = node(xk(ii));   % values of node polynomial
end
% Find position without sign changes in D*node. 
pos = find(abs(sum(sign(diag(nodevec)*Dvals))) == m+n+2 & ...
                                        sum(abs(Dvals))>1e-4);  
    
if isempty(pos) % still no solution, give up
    if ( dialogFlag && ~silentFlag )
        disp('Trial interpolant too far from optimal...')
    end
    interpSuccess = 0; 
    p = []; q = []; rh = []; h = 1e-19; wD = []; wN = [];
    return
elseif ( length(pos) > 1 ) % more than one solution with no sign changes...
    [~,ix] = min(abs(hpre)-diag(abs(d(pos,pos))));
    pos = pos(ix);
end

end    

h = -d(pos, pos);                 % levelled reference error.

% coefficients for barycentric representations
if ( m <= n )
    wD = vt(m+2:end,pos);
else
    wD = Qmn*vt(m+2:end,pos);    
end
if ( m >= n )
    wN = vt(1:m+1,pos);
else
    wN = Qmn*vt(1:m+1,pos);    
end

D = @(x) 0; N = @(x) 0;    % form function handle rh = N/D 
for ii = 1:length(xsupport)
   D = @(x) D(x) + wD(ii)./(x-xsupport(ii));
   N = @(x) N(x) + wN(ii)./(x-xsupport(ii));   
end
D = @(x)-D(x); % flip back sign

rh = @(zz) reval(zz, xsupport, N, D, wN, wD);

interpSuccess = 1; % declare success

% Form chebfuns of p and q (note: could be numerically unstable, but
% provided for convenience).
% Find values of node polynomial at Chebyshev points
if dialogFlag
    x = chebpts(m+n+1,f.domain([1,end]));
    nodex = zeros(length(x),1);
    for ii = 1:length(x)    
        nodex(ii) = node(x(ii)); 
    end 
    qvals = nodex.*feval(D,x);  % Values of p and q at Chebyshev points
    pvals = nodex.*feval(N,x);
    % If certain Chebyshev points map to support points, inf values will
    % get propagated, so we need to handle them separately
    for ii = 1:length(xsupport)
        for jj = 1:length(x)
            if x(jj) == xsupport(ii)
                nodei = 1.0;
                for kk = 1:length(xsupport)
                    if ~(kk == ii)
                        nodei = nodei * (x(jj) - xsupport(kk)); 
                    end
                end
                qvals(jj) = -nodei*wD(ii);
                pvals(jj) = nodei*wN(ii);
            end
        end
    end
 
    p = chebfun(pvals,f.domain([1,end]));
    q = chebfun(qvals,f.domain([1,end]));
    p = simplify(p); q = simplify(q);
else
    p = [];
    q = [];
end

end

function r = reval(zz, xsupport, N, D, wN, wD)
zv = zz(:);
r = N(zv)./D(zv);
ii = find(isnan(r));
for jj = 1:length(ii)
    if ( isnan(zv(ii(jj))) || ~any(zv(ii(jj)) == xsupport) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        pos = zv(ii(jj)) == xsupport; 
        r(ii(jj)) = -wN(pos)./wD(pos);
    end    
end
r = reshape(r, size(zz));
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


function [xk, norme, err_handle, flag] = exchange(xk, h, method, f, ...
                                                  fHandle, p, rh, Npts, n)
%EXCHANGE   Modify an equioscillation reference using the Remez algorithm.
%   EXCHANGE(XK, H, METHOD, F, P, RH, NPTS, N) performs one step of the
%   Remez algorithm for the best rational approximation of the CHEBFUN F
%   of the target function according to the first method (METHOD = 1),
%   i.e., exchanges only one point, or the second method (METHOD = 2),
%   i.e., exchanges all the reference points. XK is a column vector with
%   the reference, H is the levelled error, P is the numerator, and RH is a
%   function handle, NPTS is the required number of alternation points,
%   and N is the denominator degree.
%
%   [XK, NORME, E_HANDLE, FLAG] = EXCHANGE(...) returns the modified
%   reference XK, the supremum norm of the error NORME (included as an
%   output argument, since it is readily computed in EXCHANGE and is used
%   later in MINIMAX), a function handle E_HANDLE for the error, and a FLAG
%   indicating whether there were at least N+2 alternating extrema of the
%   error to form the next reference (FLAG = 1) or not (FLAG = 0).

% Compute extrema of the error.
if(n == 0) % polynomial case
    % Function handle output for evaluating the error.
    rh = @(x) feval(p,x);
    rr = findExtrema(f, fHandle, rh, xk);
    err_handle = @(x) fHandle(x) - feval(p, x);
else       % Rational case.
    rr = findExtrema(f, fHandle, rh, xk);
    err_handle = @(x) fHandle(x) - rh(x);
end

% Select exchange method.
if ( method == 1 )                             % One-point exchange.
    [~, pos] = max(abs(feval(err_handle, rr)));
    pos = pos(1);
else                                           % Full exchange.
    pos = find(abs(err_handle(rr)) >= abs(h)); % Values above levelled
                                               % error.
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
        % Given adjacent points with the same sign, keep one with largest
        % value.
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

function rts = findExtrema(f, fHandle, rh,xk)
% Finds all the local maxima and minima of f-rh.
% xk is the present reference
% rh is a handle to p/q

err_handle = @(x) fHandle(x) - rh(x);

doms = unique([f.domain'; xk]).';
doms = sort(doms,'ascend');

% Initial trial
if ( isempty(xk) )
sample_points = linspace(f.domain(1),f.domain(end),5000);
scale_of_error = norm(err_handle(sample_points),inf);
relTol =  1e-15 * (vscale(f)/scale_of_error);   
    warning off
    ek = chebfun(@(x) err_handle(x), f.domain, 'eps', relTol, ...
                                               'splitting', 'on');
    warning on
    rts = roots(diff(ek), 'nobreaks');
else
    nn = 2^3; % sampling pts in each subinterval (try low number first)
    mid = (doms(1:end-1)+doms(2:end))/2; % midpoints
    rad = (doms(2:end)-doms(1:end-1))/2; % radius
    xx = ones(nn+1,1)*mid+cos(pi*((nn:-1:0).')/nn)*rad; % sample matrix
    valerr = feval(f,xx)-rh(xx);
    rts = zeros(5*length(xk),1);
    pos = 1;
    for k = 1:length(doms)-1                
       rnow = rootsdiff(valerr(:,k),[doms(k) doms(k+1)],err_handle);
       rts(pos:pos+length(rnow)-1) = rnow; % update reference points
       pos = pos+length(rnow);
    end    
    rts(pos:end) = [];
end

% Append end points of the domain.
rts = unique([f.domain' ; rts]);
rts = sort(rts,'ascend');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find extrema of error function in AAA-Lawson
function xk = findReference(f,fHandle,r,m,n,z) 
    % f: function  
    % r: rational approximant 
    % m,n: (m,n) is the type of our rational approximant
    % z: barycentric support points 
    % OUTPUT: reference points xk
    
    % Find extrema points as usual.
    xk = findExtrema(f,fHandle, r, sort(z,'ascend'));
    
    % Deal with length(xk) not equal to the desired m+n+2.
    if length(xk) > m+n+2 % Reduce reference pts becuse too many found. 
   
        xkdiff = diff(xk);                        
        [~,ix] = sort(xkdiff,'descend'); % Take those with largest gaps.
        xk = [xk(1);xk(1+ix(1:m+n+1))];
        xk = sort(xk,'ascend');        
    elseif length(xk) < m+n+2 % Increase # of reference pts
                              % if too few found. 

        xkdiff = diff(xk);                        
        add = m+n+2-length(xk);  % We need to add this many reference
                                 % points. 
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
err = norm(err_handle(xk),'inf');
ylim(2*[-err,err])
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



function [r, pol, res, zer, z, Z, f, w, wf, errvec, p, q] = ...
                                                aaamn_lawson(F, varargin)
%AAAMN_Lawson   near-best rational approximation of F. 
% 
% R = aaamn_lawson(F) computes a rational aproximant of (default) type
% (10,10) on the default interval [-1,1]
%
% [R, POL, RES, ZER] = aaamn_lawson(...) outputs the poles, residues and
% zeros of R. The poles, zeros will approximate those of F (not well if
% R-F is not small)
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
% [...] = aaamn_lawson(F,m,n,'dom',[-1,2]) specifies the domain (this has
% no effect if Z is specified)
%
% [...] = aaamn_lawson(F,m,n,'tol',1e-5) specifies the Lawson iterate
% relative stopping criterion (here to 1e-5)
%
% [...] = aaamn_lawson(F,m,n,'iter',10) limits the Lawson maximum
% iterations to 10. 
% 
% 
% The algorithm first finds a AAA rational approximant to F \approx r(Z),
% then attempts to refine the approximant by a Lawson process, i.e. an
% iterative reweighting. 
% 
% This code is designed for computing good reference points for the
% rational minimax code to follow, but can be used independently for
% constructing a rational approximation r that can be much closer than AAA
% to the best rational approximant. 
%
% Input:  Z = vector of sample points
%         F = vector of data values, or a function handle
%         m, n: max type is (m,n), set to 10 if omitted
%         tol = relative tolerance tol, default: 1e-13 
%         Lawsoniter: max. iteration number of Lawson updates (default 10)
%         doplot: 1 to plot error curve history (default 0)
%     R = AAAMN_Lawson(F, Z, m, n, NAME, VALUE) sets the following
%     parameters:
%           - 'tol', TOL: relative tolerance for Lawson
%           iteration (default 1e-5)
%           - 'iter', IT: maximal number of Lawson iterations
%           (default MMAX = max([5 min([20, m, n])])).
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
%    [r, pol, res, zer, z, f, w, wf, errvec, p, q] = aaamn_lawson(@abs,...
%                                        10,10,'plot','on','dom',[-1 2])
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

M = length(Z);                             % number of sample points
mmax = m+1; nmax = n+1;                    % for coding convenience
if ( (nargin < 6) || isempty(Lawsoniter) ) % number of Lawson updates
    Lawsoniter = max([5 min([20,mmax,nmax])]); 
end 
if ~isfloat(F), F = feval(F,Z); end        % convert function handle to
                                           % vector
 Z = Z(:); F = F(:);                       % work with column vectors
 SF = spdiags(F,0,M,M);                    % left scaling matrix
 J = 1:M;                                  % indices that are not support
                                           % pts
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
    
        if ( mn > min(nmax,mmax) ) % nondiagonal case, find projection
                                   % subspace 
             if mmax < nmax
             q = f(:);
             else
             q = ones(length(z),1);
             end
             Q = orthspace(z,mn-min(mmax,nmax),q);  % projection subspace 
             [~,~,V] = svd(A(J,:)*Q,0);             % SVD on projected
                                                    % subspace
             w = Q*V(:,end);
        else             
             [~,~,V] = svd(A(J,:),0);               % SVD, no projection
                                                    % needed
             w = V(:,mn);                           % weight vector             
        end     
  wf = w.*f;
  N = C*(w.*f); D = C*w;                  % numerator and denominator
  R = F; R(J) = N(J)./D(J);               % rational approximation
  err = norm(F-R,inf);
  errvec = [errvec; err];                 % max error at sample points
  if ( err < tol*norm(F,inf) ), break, end    % stop if converged
end
    r = @(zz) feval(@rrint,zz,z,w,f);            % AAA approximant as
                                                 % function handle
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
                Q = orthspace(z,mn-min(mmax,nmax),q);     % projection
                                                          % subspace                              
                A =[SF*C -C*Q];                   
              end
          else
            A =[SF*C -C];   % diagonal case
          end
          
          rate = 1;         % default Lawson rate, will shrink if not
                            % converging
          nrmincreased = 0; % initialization    
	      for it = 1:Lawsoniter
              weiold = wei; 
              wei = wei .* power(abs(F(J)-R(J)),rate); % update Lawson
                                                       % weights
              wei = wei/sum(wei);                      % normalize 
              if norm(weiold-wei)/norm(wei)< tolLawson % declare Lawson
                                                       % converged
                  break
              end
              % diagonal weight matrix
              D = spdiags(sqrt(wei),0,length(wei),length(wei));

              [~,~,V] = svd(D*A(J,:),0);  % weighted least-squares via SVD
              
          if ( mn > min(nmax,mmax) )      % deal with nondiagonal case
              if ( mn > nmax )
                w = Q*V(1:nmax,end); wf = V(nmax+1:end,end);            
              else
                w = V(1:nmax,end); wf = Q*V(nmax+1:end,end); 
              end
          else
            w = V(1:mn,end); wf = V(mn+1:2*mn,end);            
          end                      
            f = wf./w;                         % for compatibility with
                                               % interpolatory-aaa                        
            N = C*wf; D = C*w;                 % numerator and denominator               
            R = F; R(J) = N(J)./D(J);          % rational approximation
            err = norm(F-R,inf);            
            errvec = [errvec; err];            % max error at sample points                            
            if ( err < nrmbest )    % adopt best so far
                nrmbest = norm(F-R,'inf'); 
                weibest = wei;                      % store best weight
                r = @(zz) feval(@rrab,zz,z,w,wf,f); % AAA approximant as
                                                    % function handle
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

            if doplot  % plot error functions
                       % (hopefully near-equioscillating)
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
        error('CHEBFUN:aaamn_lawson:nColF', ...
            'Input chebfun must have one column.')
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
     warning(['CHEBFUN:aaamn_lawson: type (m,n) not specified,', ...
         ' default to (10,10)'])
     m = 10; n = 10; 
end

% Set defaults for other parameters:
tolLawson = 1e-5;                       % Relative tolerance for
                                        % Lawson update.
tol = 1e-15;                            % AAA tolerance
Lawsoniter = max([5 min([20, m, n])]);  % Maximum number of terms.
doplot = 0;                             % Don't plot intermediate functions
                                        % unless specified

% Check if parameters have been provided:
while ( ~isempty(varargin) )
    if ( strncmpi(varargin{1}, 'tol', 3) ||  ...
            strncmpi(varargin{1}, 'tolLawson', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tolLawson = varargin{2};   % Lawson tolerance
        end
        varargin([1, 2]) = [];

    elseif ( strncmpi(varargin{1}, 'iter', 4) || ...
            strncmpi(varargin{1}, 'maxit', 5))
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            Lawsoniter = varargin{2};  % maximum Lawson iterations
        else
        warning(['CHEBFUN:aaamn_lawson:iter unspecified,', ...
            ' use default itermax ', num2str(Lawsoniter)])
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
                    ['Given domain does not match the domain ', ...
                    'of the chebfun.\n', 'Results may be inaccurate.'])
            end
        end
        
    elseif strncmpi(varargin{1}, 'plot', 4)  % plot error functions
        if isfloat(varargin{2})
            doplot = varargin{2};
        elseif ( strncmpi(varargin{2}, 'true', 4) || ...
                strncmpi(varargin{2}, 'on', 2) )
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
        error('CHEBFUN:aaamn_lawson:UnknownF', ...
            'Input for F not recognized.')
    end
    
else
    % Z was not given.  Set flag that Z needs to be determined.
    % Also set Z and M since they are needed as output.
    % in AAA this is done adaptively. This can be done with Lawson, but
    % probably safe to take as many points as reasonably possible here. 
    Z = linspace(dom(1), dom(end), 4000).';    
end

end % End of PARSEINPUT().


% generate function handle, interpolatory mode
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
function Q = orthspace(z,dim,q)    % orthonormal projection space for (m,n)
if ( dim == 0 ), Q = eye(length(z)); end 
if ( nargin < 3 ), q = ones(length(z),1); end
                Q = q/norm(q);
                for ii = 2:dim               % orthogonal complement via
                                             % Lanczos-type process
                Qtmp = diag(z)*Q(:,end);
                Qtmp = Qtmp - Q*(Q'*Qtmp);   % orthogonalize
                Qtmp = Qtmp - Q*(Q'*Qtmp);   % orthogonalize again (CGS2)                
                Qtmp = Qtmp/norm(Qtmp);      % normalize
                Q = [Q Qtmp]; 
                end
            [Q,~] = qr(Q); Q = conj(Q(:,dim+1:end)); % desired null space
end

function r = rootsdiff(vals,dom,err_handle) % vals is either a function or
                                            % values at chebpts
% returns the roots of diff(vals) via ChebyshevU-colleague
if length(vals)<=1
if nargin<3, n = 2^5; end    
xx = chebpts(n+1,dom);    
vals = vals(xx);
else 
n = size(vals,1)-1;
end

tol = 1e-3; % no need for high tolerance
c = 1;   % initialize
while ( (abs(c(end)/c(1))>tol) && (n<=2^6) )% sample until happy
cc = fft([vals(end:-1:1);vals(2:end-1)])/n;
cc(1) = cc(1)/2;
c = real(cc(1:n+1)); % coeffs of f=err_handle in T
cU = c(2:end).*(1:length(c)-1).'; % coeffs of df in U
% simplify; no need to get full accuracy.
% Then reorder to highest coeffs first. 
len = max( find((abs(cU)/norm(cU)>1e-14)) ); cU = flipud(cU(1:len)); 
%if ( length(cU)<=1 || norm(cU./max(abs(cU)))<1e-14 ), r = []; return; end
% constant function
if ( length(cU)<=1 ), r = []; return; end
if abs(c(end)/c(1))>tol  % resample at finer grid
    n = 2*n;
    vals = feval(err_handle,(dom(1)+dom(end))/2 + ...
        cos(pi*((n:-1:0).')/n)*(dom(end)-dom(1))/2);
end
end

if length(cU)<=1, r = []; return; end % constant function

if (length(cU)==2)
    % degree 1 polynomial: just take the root
    ei = -cU(2)/(2*cU(1));
    % remove root if outside domain
    ei = ei(abs(ei)<=1+1e-7);
else
    % now construct colleague matrix for ChebyshevU
    oh = ones(len-2,1)/2;
    C = diag(oh,1) + diag(oh,-1);
    cU = -cU(2:end)/cU(1)/2;cU(2) = cU(2)+.5;
    C(1,:) = cU.';
    ei = eig(C);
    % remove irrelevant roots
    ei = real(ei(abs(imag(ei))<1e-5 & abs(ei)<=1+1e-7)); 
end
r = (dom(1)+dom(2))/2 + ei*(dom(2)-dom(1))/2; % map back to the subinterval
end

% Compute the poles zeros of the barycentric approximant with
% weights alpha and beta.
function [zer, pol] = pzeros(zj, alpha, beta, rh, m, n, dom)

if ( n == 0 || isempty(beta) )
    pol = [];
else
    l = length(beta);

    % Compute poles via generalized eigenvalue problem:
    B = eye(l+1);
    B(1,1) = 0;
    E = [0 beta.'; ones(l, 1) diag(zj)];
    pol = eig(E, B);
    % Remove zeros of denominator at infinity:
    pol = pol(~isinf(pol));

    rad = 1e-5; % radius for approximating residual

    if ( l - 1 > n )% superdiagonal case, remove irrelevant poles 
    dz = rad*exp(2i*pi*(1:4)/4);
    res = rh(bsxfun(@plus, pol, dz))*dz.'/4; % residues
    ix = find( abs(res) > 1e-10 ); % pole with suff. residues
    pol = pol(ix); 
    zerBern = abs(pol-dom(1)) + abs(pol-dom(2)); % Bernstein ellipse radius
    [~,ix] = sort( zerBern, 'ascend'); % sort wrt radius
    pol = pol( ix(1:min(n,end)) ); % choose <=n zeros with largest residues    
    end
end

if ( m == 0 && isempty(alpha) )
    zer = [];
else
    l = length(alpha);
    
    % Compute zeros via generalized eigenvalue problem:
    B = eye(l+1);
    B(1,1) = 0;
    E = [0 alpha.'; ones(l, 1) diag(zj)];
    
    zer = eig(E,B);
    % Remove zeros of numerator at infinity:
    zer = zer(~isinf(zer));
    
    % subdiagonal case, remove irrelevant zeros:
    if ( l - 1 > m )
        rad = 1e-5; % radius for approximating derivative
        dz = rad*exp(2i*pi*(1:4)/4);
        deriv = sum(abs(bsxfun(@minus,rh(zer),rh(dz))./bsxfun(@minus,zer,dz)),2);
        deriv = deriv/4;
        ix = find( abs(deriv) > 1e-10 );
        
        zer = zer(ix);
        zerBern = abs(zer-dom(1)) + abs(zer-dom(2));
        [~,ix] = sort(zerBern,'ascend');
        zer = zer( ix(1:min(m,end)) );
    end
    
end

end % End of PZEROS().

function op = str2op(op)
    % Convert string inputs to either numeric format or function_handles.
    sop = str2num(op);
    if ( ~isempty(sop) )
        op = sop;
    else
        depVar = symvar(op);
        if ( numel(depVar) ~= 1 )
            error('CHEBFUN:CHEBFUN:str2op:indepvars', ...
             'Incorrect number of independent variables in string input.');
        end
        op = eval(['@(' depVar{:} ')', op]);
    end
end
