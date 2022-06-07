function [r, pol, res, zer, zj, fj, wj, errvec, wt] = aaatrig(F, varargin)
%AAATRIG   Trigonometric AAA and AAA-Lawson (near-minimax) real or complex
%      rational approximation.
%   R = AAATRIG(F, Z) computes an trigonometric AAA rational approximant R
%   (function handle) to data F on the set of sample points Z.  The rational
%   approximant is periodic with period 2*pi. F may be given by its values at
%   Z, or as a function handle or a chebfun.  R = AAATRIG(F, Z, 'degree', N)
%   computes the minimax approximation of degree N (i.e., rational type
%   (N,N)).
%
%   [R, POL, RES, ZER] = AAATRIG(F, Z) returns vectors of poles POL, residues
%   RES, and zeros ZER of R. The poles and zeros are repeated at intervals
%   of 2*pi.
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ] = AAATRIG(F, Z) also returns the vectors of
%   support points ZJ, approximation values FJ = r(ZJ), and weights WJ of
%   the barycentric representation of R.
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC] = AAATRIG(F, Z) also returns the
%   vector of errors ||f-r||_infty in successive iteration steps of AAATRIG.
%
%   R = AAATRIG(F,Z,FORM) computes a rational approximant of type FORM.
%   FORM can either be 'odd' (default) or 'even'.
%
%   R = AAATRIG(F, Z, NAME, VALUE) sets the following parameters:
%   - 'tol', TOL: relative tolerance (default TOL = 1e-13),
%   - 'degree', N: maximal degree (default N = 99). 
%      Output rational approximant will be at most of type (N,N). 
%      Identical to 'mmax', N+1. 
%      By default, this will turn on Lawson iteration: see next paragraph. 
%   - 'mmax', MMAX: maximal number of terms in the barycentric representation
%       (default MMAX = 100). R will be of degree MMAX-1. 
%       Identical to 'degree', MMAX-1. Also turns on Lawson iteration. 
%   - 'dom', DOM: domain (default DOM = [0, 2*pi]). No effect if Z is provided.
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
%   AAATRIG attempts to find a minimax approximant of degree N by Lawson iteration.
%   This will generally be successful only if the minimax error is well
%   above machine precision, and is more reliable for complex problems than
%   real ones.  If 'degree' and 'lawson' are both specified, then exactly
%   NLAWSON Lawson steps are taken (so NLAWSON = 0 corresponds to AAA
%   approximation with no Lawson iteration).  The final weight vector WT of
%   the Lawson iteration is available with 
%   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC, WT] = AAATRIG(F, Z).
%
%   Note that R may have fewer than N poles and zeros.  This may happen,
%   for example, if N is too large, or if F is even and N is odd, or if F is
%   odd and N is even.
%
%   One can also execute R = AAATRIG(F), with no specification of a set Z.
%   If F is a vector, this is equivalent to R = AAATRIG(F, Z) with
%   Z = LINSPACE(0, 2*pi, LENGTH(F)).  If F is a function handle or a chebfun,
%   AAATRIG attempts to resolve F on its domain, which defaults to [0,2*pi] for
%   a function handle.
%
% Examples:
%
%   f = @(x) exp(cos(x));
%   r = aaatrig(f); xx = linspace(-pi,pi); plot(xx,r(xx)-f(xx))
%
%   r = aaatrig(f,'degree',4); xx = linspace(-pi,pi); plot(xx,r(xx)-f(xx))
%
%   Z = exp(2i*pi*linspace(0,1,500)); 
%   [r,pol,res] = aaatrig(@tan,Z); disp([pol res])
%
%   X = linspace(0,2*pi,1000); F = gamma(2+sin(X));
%   subplot(1,2,1)
%   r = aaatrig(F,X,'degree',10,'lawson',0); plot(X,F-r(X)), hold on
%   r = aaatrig(F,X,'degree',10); plot(X,F-r(X)), hold off
% 
%   Z = exp(1i*linspace(0,2*pi,1000)); G = exp(1i*tan(Z));
%   subplot(1,2,2)
%   r = aaatrig(G,Z,'degree',5,'lawson',0); plot(G-r(Z)), axis equal, hold on
%   r = aaatrig(G,Z,'degree',5); plot(G-r(Z)), axis equal, hold off
%
%   References:
%   [1] Yuji Nakatsukasa, Olivier Sete, Lloyd N. Trefethen, "The AAA algorithm
%   for rational approximation", SIAM J. Sci. Comp. 40 (2018), A1494-A1522.
%
%   [2] Yuji Nakatsukasa and Lloyd N. Trefethen, "An algorithm for real and
%   complex rational minimax approximation", SIAM J. Sci. Comp. (2020).
%   
%   [3] Peter J. Baddoo, "The AAAtrig algorithm for rational approximation 
%   of periodic functions", SIAM J. Sci. Comp. (2021).
%
% See also AAA, TRIGRATINTERP, CHEBPADE, MINIMAX, PADEAPPROX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Parse inputs:
[F, Z, M, form, dom, tol, mmax, cleanup_flag, cleanup_tol, needZ, mmax_flag, nlawson] ...
    = parseInputs(F, varargin{:});

if ( needZ )
    % Z was not provided.  Try to resolve F on its domain.
    [r, pol, res, zer, zj, fj, wj, errvec] = ...
        aaatrig_autoZ(F, form, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag, nlawson);
    return
end

% Remove any infinite or NaN function values (avoid SVD failures):
toKeep = ~isinf(F);
F = F(toKeep); Z = Z(toKeep);
toKeep = ~isnan(F);
F = F(toKeep); Z = Z(toKeep);

% Find sample points at infity
infP = find(Z==+1i*Inf);
infM = find(Z==-1i*Inf);
finfP = F(infP); % value at +1i*Inf
finfM = F(infM); % value at -1i*Inf
F([infP,infM])=[]; Z([infP,infM])=[];
if strcmp(form,'even') && numel([finfP,finfM])==2 && finfP~=finfM
    error('The even representation must take the same values at +/-i*Inf.')
end

% Define basis functions
if strcmp(form,'even')
    cst = @(x) cot(x);
elseif strcmp(form,'odd')
    cst = @(x) csc(x);
end

% Project sample points onto a single period window
Z = Z - 2*pi*floor(real(Z/(2*pi)));

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
    C = [C cst((Z - Z(jj))/2)];
    % Compute weights:
    Sf = diag(fj);                      % Right scaling matrix.
    A = SF*C - C*Sf;                    % Loewner matrix.
    % Add value(s) at infinity to least-squares problem 
    Jv = J;
    if ~isempty([finfP,finfM]) 
        if strcmp(form,'odd')
            if ~isempty(finfP)
                A = [A; (finfP - fj.').*exp(-1i*zj.'/2)];
                Jv = [Jv, size(A,1)];
            end
            if ~isempty(finfM)
                A = [A; (finfM - fj.').*exp(+1i*zj.'/2)];
                Jv = [Jv, size(A,1)];
            end
        elseif strcmp(form,'even')
        A = [A; (finfP-fj.')]; 
        Jv = [Jv, size(A,1)];
        end
    end
    [~, ~, V] = svd(A(Jv,:), 0); % Reduced SVD. Includes the bottom
    wj = V(:,m);                        % weight vector = min sing vector
    % Rational approximant on Z:
    N = C*(wj.*fj);                     % Numerator
    D = C*wj;                           % Denominator
    R = F;
    R(J) = N(J)./D(J);
    
    % Error in the sample points and at infinity:
    err = F - R;
    if isempty([finfP,finfM])
    maxerr = norm(1, inf);
    else
        if strcmp(form,'odd')
            if ~isempty(finfP)
            errInfP = finfP - (fj.'.*exp(-1i*zj.'/2)*wj)./(exp(-1i*zj.'/2)*wj); % error at +1i*infinity
            err = [err; errInfP];
            end
            if ~isempty(finfM)
            errInfM = finfM - (fj.'.*exp(+1i*zj.'/2)*wj)./(exp(+1i*zj.'/2)*wj); % error at -1i*infinity
            err = [err; errInfM];
            end
        elseif strcmp(form,'even')
            errInf = finfP - fj.'*wj./sum(wj); % error at +/-1i*infinity
            err = [err; errInf];
        end
    end
    maxerr = norm(err,inf);
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
    wj = [1; -1];       
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

  if ~isempty([finfP,finfM])
    warning('Specifying the function values at infinity is not currently compatible with Lawson iteration.')
  end
   
    maxerrold = maxerrAAA;
    maxerr = maxerrold;
    nj = length(zj);
    A = [];
    for j = 1:nj                              % Cauchy/Loewner matrix
        A = [A cst((Z-zj(j))/2) F.*cst((Z-zj(j))/2)];
    end
    for j = 1:nj
        [i,~] = find(Z==zj(j));               % support pt rows are special
        A(i,:) = 0;
        A(i,2*j-1) = 2;
        A(i,2*j) = 2*F(i);
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
            denom = denom + c(2*j).*cst((Z-zj(j))/2);
            num = num - c(2*j-1).*cst((Z-zj(j))/2);
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
r = @(zz) revaltrig(zz, zj, fj, wj, form);

% Compute poles, residues and zeros:
[pol, res, zer] = prztrig(zj, fj, wj, form);

if ( cleanup_flag & nlawson == 0)       % Remove Froissart doublets
    [r, pol, res, zer, zj, fj, wj] = ...
        cleanuptrig(r, form, finfP, finfM,pol, res, zer, zj, fj, wj, Z, F, cleanup_tol);
end

end % of AAATRIG()

%% parse Inputs:

function [F, Z, M, form, dom, tol, mmax, cleanup_flag, cleanup_tol, ...
    needZ, mmax_flag, nlawson] = parseInputs(F, varargin)
% Input parsing for AAATRIG.

% Check if F is empty:
if ( isempty(F) )
    error('CHEBFUN:aaatrig:emptyF', 'No function given.')
elseif ( isa(F, 'chebfun') )
    if ( size(F, 2) ~= 1 )
        error('CHEBFUN:aaatrig:nColF', 'Input chebfun must have one column.')
    end
end

% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    % Z is given.
    Z = varargin{1};
    if ( isempty(Z) )
        error('CHEBFUN:aaatrig:emptyZ', ...
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
    dom = [0, 2*pi];
end
cleanup_flag = 1;   % Cleanup on.
mmax_flag = 0;      % Checks if mmax manually specified.
cleanup_set = 0;    % Checks if cleanup_tol manually specified.
form = 'odd';
% Check if parameters have been provided:
while ( ~isempty(varargin) )
    
    if  strncmpi(varargin{1},'even',4)
          form = 'even';
          varargin(1) = [];
    elseif strncmpi(varargin{1},'odd',3)
          form = 'odd';
          varargin(1) = [];
        
    elseif ( strncmpi(varargin{1}, 'tol', 3) )
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
                error('CHEBFUN:aaatrig:degmmaxmismatch', ' mmax must equal degree+1.')
            end            
            mmax = varargin{2}+1;
            mmax_flag = 1;
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )            
            if ( mmax_flag == 1 ) && ( mmax ~= varargin{2})                
                error('CHEBFUN:aaatrig:degmmaxmismatch', ' mmax must equal degree+1.')
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
                warning('CHEBFUN:aaatrig:dom', ...
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
        end
        varargin([1, 2]) = [];
        
    else
        error('CHEBFUN:aaatrig:UnknownArg', 'Argument unknown.')
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
            error('CHEBFUN:aaatrig:lengthFZ', ...
                'Inputs F and Z must have the same length.')
        end
    elseif ( ischar(F) )
        % F is given as a string input. Convert it to a function handle.
        F = str2op(vectorize(F));
        F = F(Z);
    else
        error('CHEBFUN:aaatrig:UnknownF', 'Input for F not recognized.')
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

end % End of PARSEINPUT().


%% Cleanup.  In June 2022 the residue size test was changed to be relative to
%            the distance to the approximation set Z.

function [r, pol, res, zer, z, f, w] = ...
    cleanuptrig(r, form, finfP, finfM, pol, res, zer, z, f, w, Z, F, cleanup_tol) 
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
    warning('CHEBFUN:aaatrig:Froissart','1 Froissart doublet');
else
    warning('CHEBFUN:aaatrig:Froissart',[int2str(ni) ' Froissart doublets']);
end
% For each spurious pole find and remove closest support point.
for j = 1:ni
    % Find the closest support point modulo the period
    np = fix(real((z-pol(ii(j)))/pi));
    azp = abs(z-(pol(ii(j)) + np*2*pi));
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

% Define basis functions
if strcmp(form,'even')
    cst = @(x) cot(x);
elseif strcmp(form,'odd')
    cst = @(x) csc(x);
end

% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(f);
C = cst(bsxfun(@minus, Z, z.')/2);      % Cauchy matrix.
A = SF*C - C*Sf;                    % Loewner matrix.

% Enforce value(s) at infinity 
if ~isempty([finfP,finfM]) 
    if strcmp(form,'odd')
        if ~isempty(finfP)
            A = [A; (finfP - f.').*exp(-1i*z.'/2)];
        end
        if ~isempty(finfM)
            A = [A; (finfM - f.').*exp(+1i*z.'/2)];
        end
    elseif strcmp(form,'even')
    A = [A; (finfP-f.')];   % Add a row to A to enforce behaviour at infinity.
    end
end

% Solve least-squares problem to obtain weights:
[~, ~, V] = svd(A, 0);
w = V(:,m);

% Build function handle and compute poles, residues and zeros:
r = @(zz) revaltrig(zz, z, f, w, form);
[pol, res, zer] = prztrig(z, f, w, form);

end % End of CLEANUPTRIG().


%% Automated choice of sample set

function [r, pol, res, zer, zj, fj, wj, errvec] = ...
    aaatrig_autoZ(F, form, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag, nlawson)
%

% Flag if function has been resolved:
isResolved = 0;

% Main loop:
for n = 5:14
    % Sample points:
    % Next line enables us to do pretty well near poles
    Z = linspace(dom(1)+1.37e-7*diff(dom), dom(2), 2 + 2^n).'; Z(end) = [];
    [r, pol, res, zer, zj, fj, wj, errvec] = aaatrig(F, Z, form, 'tol', tol, ...
        'mmax', mmax, 'cleanup', cleanup_flag, 'cleanuptol', cleanup_tol, 'lawson', nlawson);
    
    % Test if rational approximant is accurate:
    reltol = tol * norm(F(Z), inf);
    
    % On Z(n):
    err(1,1) = norm(F(Z) - r(Z), inf);
    Zrefined = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), ...
        1 + round(1.5 * (1 + 2^(n+1)))).';
    Zrefined(end) = [];
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
    warning('CHEBFUN:aaatrig:notResolved', ...
        'Function not resolved using %d pts.', length(Z))
end

end % End of AAATRIG_AUTOZ().

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
