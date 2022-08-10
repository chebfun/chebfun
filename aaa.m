function [r, pol, res, zer, zj, fj, wj, errvec, wt] = aaa(F, varargin)
%AAA   AAA and AAA-Lawson (near-minimax) real or complex rational approximation.
%   R = AAA(F, Z) computes the AAA rational approximant R (function handle) to
%   data F on the set of sample points Z.  F may be given by its values at Z,
%   or as a function handle or a chebfun.  R = AAA(F, Z, 'degree', N) computes
%   the minimax approximation of degree N (i.e., rational type (N,N)).
%
%   [R, POL, RES, ZER] = AAA(F, Z) returns vectors of poles POL, residues RES,
%   and zeros ZER of R.
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ] = AAA(F, Z) also returns the vectors
%   of support points ZJ, approximation values FJ = r(ZJ), and weights WJ 
%   of the barycentric representation of R. 
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC] = AAA(F, Z) also returns the
%   vector of errors ||f-r||_infty in successive iterative steps of AAA.
%   Note that the rational degrees are not 1,2,...,length(ERRVEC) but
%   0,1,...,length(ERRVEC)-1.
%
%   R = AAA(F, Z, NAME, VALUE) sets the following parameters:
%   - 'tol', TOL: relative tolerance (default TOL = 1e-13),
%   - 'degree', N: maximal degree (default N = 99). 
%      Output rational approximant will be at most of type (N,N). 
%      Identical to 'mmax', N+1. 
%      By default, this will turn on Lawson iteration: see next paragraph. 
%   - 'mmax', MMAX: maximal number of terms in the barycentric representation
%       (default MMAX = 100). R will be of degree MMAX-1. 
%       Identical to 'degree', MMAX-1. Also turns on Lawson iteration. 
%   - 'dom', DOM: domain (default DOM = [-1, 1]). No effect if Z is provided.
%   - 'cleanup', 'off' or 0: turns off automatic removal of numerical Froissart
%       doublets
%   - 'cleanuptol', CLEANUPTOL: cleanup tolerance (default CLEANUPTOL = TOL).
%       Poles with residues less than this number times the geometric mean size
%       of F times the minimum distance to Z are deemed spurious by the cleanup
%       procedure. If TOL = 0, then CLEANUPTOL defaults to 1e-13.
%   - 'lawson', NLAWSON: take NLAWSON iteratively reweighted least-squares steps
%       to bring approximation closer to minimax; specifying NLAWSON = 0 
%       ensures there is no Lawson iteration.  See next paragraph.
%
%   If 'degree' or equivalently 'mmax' is specified and 'lawson' is not, then
%   AAA attempts to find a minimax approximant of degree N by Lawson iteration.
%   This will generally be successful only if the minimax error is well
%   above machine precision, and is more reliable for complex problems than
%   real ones.  If 'degree' and 'lawson' are both specified, then exactly
%   NLAWSON Lawson steps are taken (so NLAWSON = 0 corresponds to AAA
%   approximation with no Lawson iteration).  The final weight vector WT of
%   the Lawson iteration is available with 
%   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC, WT] = AAA(F, Z).
%
%   Note that R may have fewer than N poles and zeros.  This may happen,
%   for example, if N is too large, or if F is even and N is odd, or if F is
%   odd and N is even.
%
%   One can also execute R = AAA(F), with no specification of a set Z.
%   If F is a vector, this is equivalent to R = AAA(F, Z) with
%   Z = LINSPACE(-1, 1, LENGTH(F)).  If F is a function handle or a chebfun,
%   AAA attempts to resolve F on its domain, which defaults to [-1,1] for
%   a function handle.
%
% Examples:
%   r = aaa(@exp); xx = linspace(-1,1); plot(xx,r(xx)-exp(xx))
%
%   r = aaa(@exp,'degree',4); xx = linspace(-1,1); plot(xx,r(xx)-exp(xx))
%
%   Z = exp(2i*pi*linspace(0,1,500)); 
%   [r,pol,res] = aaa(@tan,Z); disp([pol res])
%
%   X = linspace(-1,1,1000); F = tanh(20*X);
%   subplot(1,2,1)
%   r = aaa(F,X,'degree',15,'lawson',0); plot(X,F-r(X)), hold on
%   r = aaa(F,X,'degree',15); plot(X,F-r(X)), hold off
% 
%   Z = exp(1i*pi*linspace(-1,1,1000)); G = exp(Z);
%   subplot(1,2,2)
%   r = aaa(G,Z,'degree',3,'lawson',0); plot(G-r(Z)), axis equal, hold on
%   r = aaa(G,Z,'degree',3); plot(G-r(Z)), axis equal, hold off
%
%   References:
%   [1] Yuji Nakatsukasa, Olivier Sete, Lloyd N. Trefethen, "The AAA algorithm
%   for rational approximation", SIAM J. Sci. Comp. 40 (2018), A1494-A1522.
%
%   [2] Yuji Nakatsukasa and Lloyd N. Trefethen, An algorithm for real and
%   complex rational minimax approximation, SIAM J. Sci. Comp. 42 (2020),
%   A3157-A3179.
%
% See also AAATRIG, CF, CHEBPADE, MINIMAX, PADEAPPROX, RATINTERP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Parse inputs:
[F, Z, M, dom, tol, mmax, cleanup_flag, cleanup_tol, needZ, mmax_flag, nlawson] ...
    = parseInputs(F, varargin{:});

if ( needZ )
    % Z was not provided.  Try to resolve F on its domain.
    [r, pol, res, zer, zj, fj, wj, errvec] = ...
        aaa_autoZ(F, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag, nlawson);
    return
end

% Remove any infinite or NaN function values (avoid SVD failures):
toKeep = ~isinf(F);
F = F(toKeep); Z = Z(toKeep);
toKeep = ~isnan(F);
F = F(toKeep); Z = Z(toKeep);

% Remove repeated elements of Z and corresponding elements of F:
[Z, uni] = unique(Z,'stable'); F = F(uni);

M = length(Z);

% Relative tolerance:
reltol = tol * norm(F, inf);

% Left scaling matrix:
SF = spdiags(F, 0, M, M);

% Initialization for AAA iteration:
J = 1:M;
zj = [];
fj = [];
C = [];
errvec = [];
R = mean(F);

% AAA iteration:
for m = 1:mmax
    % Select next support point where error is largest:
    [~, jj] = max(abs(F - R));          % Select next support point.
    zj = [zj; Z(jj)];                   % Update support points.
    fj = [fj; F(jj)];                   % Update data values.
    J(J == jj) = [];                    % Update index vector.
    C = [C 1./(Z - Z(jj))];             % Next column of Cauchy matrix.
    
    % Compute weights:
    Sf = diag(fj);                      % Right scaling matrix.
    A = SF*C - C*Sf;                    % Loewner matrix.
    [~, ~, V] = svd(A(J,:), 0);         % Reduced SVD.
    wj = V(:,m);                        % weight vector = min sing vector
    
    % Rational approximant on Z:
    N = C*(wj.*fj);                     % Numerator
    D = C*wj;                           % Denominator
    R = F;
    R(J) = N(J)./D(J);
    
    % Error in the sample points:
    maxerr = norm(F - R, inf);
    errvec = [errvec; maxerr];
    
    % Check if converged:
    if ( maxerr <= reltol )
        break
    end
end
maxerrAAA = maxerr;                     % error at end of AAA 

% When M == 2, one weight is zero and r is constant.
% To obtain a good approximation, interpolate in both sample points.
if ( M == 2 )
    zj = Z;
    fj = F;
    wj = [1; -1];       % Only pole at infinity.
    wj = wj/norm(wj);   % Impose norm(w) = 1 for consistency.
    errvec(2) = 0;
    maxerrAAA = 0;
end

% We now enter Lawson iteration: barycentric IRLS = iteratively reweighted
% least-squares if 'lawson' is specified with NLAWSON > 0 or 'mmax' is
% specified and 'lawson' is not.  In the latter case the number of steps
% is chosen adaptively.  Note that the Lawson iteration is unlikely to be
% successful when the errors are close to machine precision.

wj0 = wj; fj0 = fj;     % Save parameters in case Lawson fails
wt = NaN(M,1); wt_new = ones(M,1);
if ( nlawson > 0 )      % Lawson iteration

    maxerrold = maxerrAAA;
    maxerr = maxerrold;
    nj = length(zj);
    A = [];
    for j = 1:nj                              % Cauchy/Loewner matrix
        A = [A 1./(Z-zj(j)) F./(Z-zj(j))];
    end
    for j = 1:nj
        [i,~] = find(Z==zj(j));               % support pt rows are special
        A(i,:) = 0;
        A(i,2*j-1) = 1;
        A(i,2*j) = F(i);
    end
    stepno = 0;
    while ( (nlawson < inf) & (stepno < nlawson) ) |...
          ( (nlawson == inf) & (stepno < 20) ) |...
          ( (nlawson == inf) & (maxerr/maxerrold < .999) & (stepno < 1000) ) 
        stepno = stepno + 1;
        wt = wt_new;
        W = spdiags(sqrt(wt),0,M,M);
        [U,S,V] = svd(W*A,0);
        c = V(:,end);
        denom = zeros(M,1); num = zeros(M,1);
        for j = 1:nj
            denom = denom + c(2*j)./(Z-zj(j));
            num = num - c(2*j-1)./(Z-zj(j));
        end
        R = num./denom;
        for j = 1:nj
            [i,~] = find(Z==zj(j));           % support pt rows are special
            R(i) = -c(2*j-1)/c(2*j);
        end
        err = F - R; abserr = abs(err);
        wt_new = wt.*abserr; wt_new = wt_new/norm(wt_new,inf);
        maxerrold = maxerr;
        maxerr = max(abserr);
    end
    wj = c(2:2:end);
    fj = -c(1:2:end)./wj;
    % If Lawson has not reduced the error, return to pre-Lawson values.
    if (maxerr > maxerrAAA) & (nlawson == Inf)
        wj = wj0; fj = fj0; 
    end
end

% Remove support points with zero weight:
I = find(wj == 0);
zj(I) = [];
wj(I) = [];
fj(I) = [];

% Construct function handle:
r = @(zz) reval(zz, zj, fj, wj);

% Compute poles, residues and zeros:
[pol, res, zer] = prz(zj, fj, wj);

if ( cleanup_flag == 1 && nlawson == 0)       % Remove Froissart doublets
    [r, pol, res, zer, zj, fj, wj] = ...
        cleanup(r, pol, res, zer, zj, fj, wj, Z, F, cleanup_tol);
elseif ( cleanup_flag == 2 && nlawson == 0)   % Alternative cleanup.  For the
                                              % moment this is an undocumented
                                              % feature, pending further
                                              % investigation.
    a.zj = zj; a.fj = fj; a.wj = wj;
    a.Z = Z; a.F = F;
    a.cleanup_tol = max(cleanup_tol, eps);
    c = cleanup2(a);
    zj = c.zj; fj = c.fj; wj = c.wj;
    r = @(zz) reval(zz, zj, fj, wj);
    [pol, res, zer] = prz(zj, fj, wj);
end

end % of AAA()


%% parse Inputs:

function [F, Z, M, dom, tol, mmax, cleanup_flag, cleanup_tol, ...
    needZ, mmax_flag, nlawson] = parseInputs(F, varargin)
% Input parsing for AAA.

% Check if F is empty:
if ( isempty(F) )
    error('CHEBFUN:aaa:emptyF', 'No function given.')
elseif ( isa(F, 'chebfun') )
    if ( size(F, 2) ~= 1 )
        error('CHEBFUN:aaa:nColF', 'Input chebfun must have one column.')
    end
end

% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    % Z is given.
    Z = varargin{1};
    if ( isempty(Z) )
        error('CHEBFUN:aaa:emptyZ', ...
            'If sample set is provided, it must be nonempty.')
    end
    varargin(1) = [];
end

% Set defaults for other parameters:
tol = 1e-13;         % Relative tolerance.
mmax = 100;          % Maximum number of terms.
cleanup_tol = 1e-13; % Cleanup tolerance.
nlawson = Inf;       % number of Lawson steps (Inf means adaptive)
% Domain:
if ( isa(F, 'chebfun') )
    dom = F.domain([1, end]);
else
    dom = [-1, 1];
end
cleanup_flag = 1;   % Cleanup on.
mmax_flag = 0;      % Checks if mmax manually specified.
cleanup_set = 0;    % Checks if cleanup_tol manually specified.
% Check if parameters have been provided:
while ( ~isempty(varargin) )
    if ( strncmpi(varargin{1}, 'tol', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tol = varargin{2};
            if ~cleanup_set & tol > 0 % If not manually set, set cleanup_tol to tol.
              cleanup_tol = tol;
            end
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'degree', 6) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            if ( mmax_flag == 1 ) && ( mmax ~= varargin{2}+1 )
                error('CHEBFUN:aaa:degmmaxmismatch', ' mmax must equal degree+1.')
            end            
            mmax = varargin{2}+1;
            mmax_flag = 1;
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )            
            if ( mmax_flag == 1 ) && ( mmax ~= varargin{2})                
                error('CHEBFUN:aaa:degmmaxmismatch', ' mmax must equal degree+1.')
            end
            mmax = varargin{2};
            mmax_flag = 1;
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'lawson', 6) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            nlawson = varargin{2};
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'dom', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 2]) )
            dom = varargin{2};
        end
        varargin([1, 2]) = [];
        if ( isa(F, 'chebfun') )
            if ( ~isequal(dom, F.domain([1, end])) )
                warning('CHEBFUN:aaa:dom', ...
                    ['Given domain does not match the domain of the chebfun.\n', ...
                    'Results may be inaccurate.'])
            end
        end
        
    elseif ( strncmpi(varargin{1}, 'cleanuptol', 10) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
          cleanup_tol = varargin{2};
          cleanup_set = 1;
        end
        varargin([1, 2]) = [];

    elseif ( strncmpi(varargin{1}, 'cleanup', 7) )
        if ( strncmpi(varargin{2}, 'off', 3) || ( varargin{2} == 0 ) )
            cleanup_flag = 0;
        elseif ( varargin{2} == 2 )     % Alternative cleanup
            cleanup_flag = 2;
        end
        varargin([1, 2]) = [];
        
    else
        error('CHEBFUN:aaa:UnknownArg', 'Argument unknown.')
    end
end

% Deal with Z and F:
if ( ~exist('Z', 'var') && isfloat(F) )
    % F is given as data values, pick same number of sample points:
    Z = linspace(dom(1), dom(2), length(F)).';
end

if ( exist('Z', 'var') )
    % Z is given:
    needZ = 0;
    
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
            error('CHEBFUN:aaa:lengthFZ', ...
                'Inputs F and Z must have the same length.')
        end
    elseif ( ischar(F) )
        % F is given as a string input. Convert it to a function handle.
        F = str2op(vectorize(F));
        F = F(Z);
    else
        error('CHEBFUN:aaa:UnknownF', 'Input for F not recognized.')
    end
    
else
    % Z was not given.  Set flag that Z needs to be determined.
    % Also set Z and M since they are needed as output.
    needZ = 1;
    Z = [];
    M = length(Z);
end

if ~mmax_flag & (nlawson == Inf)
    nlawson = 0;               
end

end % End of PARSEINPUTS().


%% Cleanup.  In June 2022 the residue size test was changed to be relative to 
%            the distance to the approximation set Z.

function [r, pol, res, zer, z, f, w] = ...
    cleanup(r, pol, res, zer, z, f, w, Z, F, cleanup_tol) 
% Remove spurious pole-zero pairs.

% Find negligible residues:
if any(F)
   geometric_mean_of_absF = exp(mean(log(abs(F(F~=0)))));
else
   geometric_mean_of_absF = 0;
end
Zdistances = NaN(size(pol));
for j = 1:length(Zdistances);
   Zdistances(j) = min(abs(pol(j)-Z));
end
ii = find(abs(res)./Zdistances < cleanup_tol * geometric_mean_of_absF);
ni = length(ii);
if ( ni == 0 )
    % Nothing to do.
    return
elseif ( ni == 1 )
    warning('CHEBFUN:aaa:Froissart','1 Froissart doublet');
else
    warning('CHEBFUN:aaa:Froissart',[int2str(ni) ' Froissart doublets']);
end

% For each spurious pole find and remove closest support point:
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
    
    % Remove support point(s):
    z(jj) = [];
    f(jj) = [];
end

% Remove support points z from sample set:
for jj = 1:length(z)
    F(Z == z(jj)) = [];
    Z(Z == z(jj)) = [];
end
m = length(z);
M = length(Z);

% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(f);
C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix.
A = SF*C - C*Sf;                    % Loewner matrix.

% Solve least-squares problem to obtain weights:
[~, ~, V] = svd(A, 0);
w = V(:,m);

% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, z, f, w);
[pol, res, zer] = prz(z, f, w);

end % End of CLEANUP().

function c = cleanup2(a)
%% Alternative cleanup procedure to remove spurious pole-zero pairs.
% This considers pole-zero distances.  Stefano Costa, August 2022.

z = a.zj; f = a.fj; w = a.wj;
[pol, res, zer] = prz(z, f, w);

cleanup_tol = a.cleanup_tol;

niter = 0;
while(true)
    niter = niter+1;
    Z = a.Z; F = a.F;
    ii = [];
    for jj = 1:length(pol)
        dz = min(abs(zer-pol(jj))); if isempty(dz), dz = 1e100; end
        dS = abs(Z-pol(jj));
        ds = min(dS);
        if any(F)
            q = 4*pi*abs(F).*dS;
            Q = mean(q);                % Arithmetic mean
        else
            Q = 0;
        end
        R = 8*cleanup_tol*Q/(4*pi);    % Equivalent residue value
        
        % Conditions to expunge poles
        % Expunge if either minimum distance is zero
        if (ds==0) || (dz==0)
            ii = [ii; jj];
        % Expunge if Z is a real interval
        elseif isreal(Z) && (abs(imag(pol(jj)))<eps) && ...
            (real(pol(jj))>=min(Z)) && (real(pol(jj))<=max(Z))
            ii = [ii; jj];
        % Expunge if Z is the unit disk
        elseif all(abs(Z)==1) && (abs(abs(pol(jj))-1)<eps)
            ii = [ii; jj];
        % Expunge if distance to closest zero is undetectable
        elseif (dz/ds<1) && (dz<max(cleanup_tol^2,eps))
            ii = [ii; jj];
        % Expunge if a nearby zero exists and residue is below the
        % equivalent value R. Two choices for real and complex F
        elseif ((dz/ds)<sqrt(cleanup_tol))
            if ( ~any(imag(F)) && (abs(real(res(jj))) < R) )
                ii = [ii; jj];
            elseif (abs(res(jj)) < R)
                ii = [ii; jj];
            end
        end
    end
    ii = unique(ii);

    ni = length(ii);
    if ( ni == 0 )
        % Nothing to do.
        break;
    elseif ( ni == 1 )
        warning('CHEBFUN:aaa:Froissart',...
            ['1 Froissart doublet, niter = ', int2str(niter)]);
    else
        warning('CHEBFUN:aaa:Froissart',...
            [int2str(ni) ' Froissart doublets, niter = ' int2str(niter)]);
    end

    % For each spurious pole find and remove closest support point:
    for j = 1:ni
        azp = abs(z-pol(ii(j)));
        jj = find(azp == min(azp),1);
        
        % Remove support point(s):
        z(jj) = [];
        f(jj) = [];
    end
    
    % Remove support points z from sample set:
    for jj = 1:length(z)
        F(Z == z(jj)) = [];
        Z(Z == z(jj)) = [];
    end
    m = length(z);
    M = length(Z);
    
    % Build Loewner matrix:
    SF = spdiags(F, 0, M, M);
    Sf = diag(f);
    C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix.
    A = SF*C - C*Sf;                    % Loewner matrix.
    
    % Solve least-squares problem to obtain weights:
    [~, ~, V] = svd(A, 0);
    w = V(:,m);
    
    % Compute poles, residues and zeros:
    [pol, res, zer] = prz(z, f, w);
end % End of while loop

c.zj = z; c.fj = f; c.wj = w;

end  % End of CLEANUP2().

%% Automated choice of sample set

function [r, pol, res, zer, zj, fj, wj, errvec] = ...
    aaa_autoZ(F, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag, nlawson)
%

% Flag if function has been resolved:
isResolved = 0;

% Main loop:
for n = 5:14
    % Sample points:
    % Next line enables us to do pretty well near poles
    Z = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), 1 + 2^n).';
    [r, pol, res, zer, zj, fj, wj, errvec] = aaa(F, Z, 'tol', tol, ...
        'mmax', mmax, 'cleanup', cleanup_flag, 'cleanuptol', cleanup_tol, 'lawson', nlawson);
    
    % Test if rational approximant is accurate:
    reltol = tol * norm(F(Z), inf);
    
    % On Z(n):
    err(1,1) = norm(F(Z) - r(Z), inf);
    
    Zrefined = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), ...
        round(1.5 * (1 + 2^(n+1)))).';
    err(2,1) = norm(F(Zrefined) - r(Zrefined), inf);
    
    if ( all(err < reltol) )
        % Final check that the function is resolved, inspired by sampleTest().
        % Pseudo random sample points in [-1, 1]:
        xeval = [-0.357998918959666; 0.036785641195074];
        % Scale to dom:
        xeval = (dom(2) - dom(1))/2 * xeval + (dom(2) + dom(1))/2;
        
        if ( norm(F(xeval) - r(xeval), inf) < reltol )
            isResolved = 1;
            break
        end
    end
end

if ( ( isResolved == 0 ) && ~mmax_flag )
    warning('CHEBFUN:aaa:notResolved', ...
        'Function not resolved using %d pts.', length(Z))
end

end % End of AAA_AUTOZ().

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
end % End of STR2OP().
