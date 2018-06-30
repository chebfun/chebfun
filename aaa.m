function [r, pol, res, zer, zj, fj, wj, errvec] = aaa(F, varargin)
%AAA   Computes a AAA rational approximation.
%   R = AAA(F, Z) computes the AAA rational approximant R (function handle) to
%   data F on the set of sample points Z.  F may be given by its values at Z,
%   or as a function handle or a chebfun.
%
%   [R, POL, RES, ZER] = AAA(F, Z) returns vectors of poles POL,
%   residues RES, and zeros ZER of R.
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ] = AAA(F, Z) also returns the vectors
%   of support points ZJ, function values FJ, and weights WJ of the
%   barycentric representation of R.
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC] = AAA(F, Z) also returns the
%   vector of errors ||f-r||_infty in successive iteration steps of AAA.
%
%   R = AAA(F, Z, NAME, VALUE) sets the following parameters:
%   - 'tol', TOL: relative tolerance (default TOL = 1e-13),
%   - 'mmax', MMAX: maximal number of terms in the barycentric representation
%       (default MMAX = 100).
%   - 'dom', DOM: domain (default DOM = [-1, 1]). No effect if Z is provided.
%   - 'cleanup', 'off' or 0: turns off automatic removal of numerical Froissart
%       doublets
%   - 'cleanuptol', CLEANUPTOL: cleanup tolerance (default CLEANUPTOL = TOL).
%       Poles with residues less than this number times the maximium absolute
%       component of F are deemed spurious by the cleanup procedure. If TOL = 0,
%       then CLEANUPTOL defaults to 1e-13.
%
%   One can also execute R = AAA(F), with no specification of a set Z.
%   This is equivalent to defining Z = LINSPACE(DOM(1), DOM(2), LENGTH(F)) if F
%   is a vector (by default DOM = [-1, 1]).
%   If F is a function handle or a chebfun, AAA attempts to resolve F on its
%   domain DOM.  By default, DOM = [-1, 1] for a function handle, and
%   DOM = F.DOMAIN([1, END]) for a chebfun.
%
% Examples:
%   r = aaa(@exp); xx = linspace(-1,1); plot(xx,r(xx)-exp(xx))
%
%   Z = exp(2i*pi*linspace(0,1,500)); 
%   [r,pol,res] = aaa(@tan,Z); disp([pol res])
%
%   Reference:
%   [1] Yuji Nakatsukasa, Olivier Sete, Lloyd N. Trefethen, "The AAA algorithm
%   for rational approximation", SIAM J. Sci. Comp. 40 (2018), A1494-A1522.
%
% See also CF, CHEBPADE, MINIMAX, PADEAPPROX, RATINTERP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Parse inputs:
[F, Z, M, dom, tol, mmax, cleanup_flag, cleanup_tol, needZ, mmax_flag, nlawson] ...
    = parseInputs(F, varargin{:});

if ( needZ )
    % Z was not provided.  Try to resolve F on its domain.
    [r, pol, res, zer, zj, fj, wj, errvec] = ...
        aaa_autoZ(F, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag);
    return
end

% Remove any infinite or NaN function values (avoid SVD failures):
toKeep = ~isinf(F);
F = F(toKeep); Z = Z(toKeep);
toKeep = ~isnan(F);
F = F(toKeep); Z = Z(toKeep);

% Remove repeated elements of Z and corresponding elements of F:
[Z, uni] = unique(Z); F = F(uni);

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
    err = norm(F - R, inf);
    errvec = [errvec; err];
    
    % Check if converged:
    if ( err <= reltol )
        break
    end
end

% Note: When M == 2, one weight is zero and r is constant.
% To obtain a good approximation, interpolate in both sample points.
if ( M == 2 )
    zj = Z;
    fj = F;
    wj = [1; -1];       % Only pole at infinity.
    wj = wj/norm(wj);   % Impose norm(w) = 1 for consistency.
    errvec(2) = 0;
end

% The following Lawson iteration option is, as of mid-2018,
% an undocumented feature.  The reason for introducing this is
% that we are in the midst of extensive experiments with 
% the AAA-Lawson idea, for which we want it to be conveniently
% available as an option in the aaa code.  However, until the
% experiments reach their conclusion, we don't regard the Lawson
% idea as settled and reliable enough to be publicized to users
% in the help text.  The default operation of aaa invokes no
% Lawson steps, so users will not notice the existence of this
% feature.
%
% Note that Lawson steps are generally useless for approximations
% down near machine precision.  They should be used in the
% context of a looser aaa 'tol' or a restricted aaa 'mmax'.
% 
% Examples:
%
% X = linspace(-1,1,1000); F = tanh(20*X);
% subplot(1,2,1)
% r = aaa(F,X,'mmax',16); plot(X,F-r(X)), hold on
% r = aaa(F,X,'mmax',16,'lawson',50); plot(X,F-r(X)), hold off
% 
% Z = exp(1i*pi*linspace(-1,1,1000)); G = exp(Z);
% subplot(1,2,2)
% r = aaa(G,Z,'tol',1e-6); plot(G-r(Z)), axis equal, hold on
% r = aaa(G,Z,'tol',1e-6,'lawson',10); plot(G-r(Z)), axis equal, hold off
% 
%         - Nick Trefethen and Abi Gopal, 29 June 2018.

if ( nlawson > 0 )      % nlawson steps of Lawson iteration

  nj = length(zj); Z2 = Z; F2 = F; np2 = M - nj;
  for j = 1:nj
    [i,~] = find(Z2==zj(j));
    Z2(i) = []; F2(i) = [];
  end
  A = []; for j = 1:nj, A = [A 1./(Z2-zj(j)) F2./(Z2-zj(j))]; end
  wt = ones(np2,1);
  for n = 0:nlawson
    W = spdiags(sqrt(wt),0,np2,np2); [U,S,V] = svd(W*A,0); c = V(:,end);
    denom = zeros(np2,1); num = zeros(np2,1);
    for j = 1:nj
      denom = denom + c(2*j)./(Z2-zj(j));
      num = num - c(2*j-1)./(Z2-zj(j));
    end
    err = F2 - num./denom; abserr = abs(err);
    wt = wt.*abserr; wt = wt/norm(wt,inf);
  end
  wj = c(2:2:end);
  fj = -c(1:2:end)./wj;

end

% Remove support points with zero weight:
I = find(wj == 0);
zj(I) = [];
wj(I) = [];
fj(I) = [];

% Construct function handle:
r = @(zz) reval(zz, zj, fj, wj);

% Compute poles, residues and zeros:
[pol, res, zer] = prz(r, zj, fj, wj);

if ( cleanup_flag )
    % Remove Froissart doublets:
    [r, pol, res, zer, zj, fj, wj] = ...
        cleanup(r, pol, res, zer, zj, fj, wj, Z, F, cleanup_tol);
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
nlawson = 0;         % number of Lawson steps
% Domain:
if ( isa(F, 'chebfun') )
    dom = F.domain([1, end]);
else
    dom = [-1, 1];
end
cleanup_flag = 1;   % Cleanup on.
mmax_flag = 0;
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
        
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
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

end % End of PARSEINPUT().


%% Cleanup

function [r, pol, res, zer, z, f, w] = ...
    cleanup(r, pol, res, zer, z, f, w, Z, F, cleanup_tol) 
% Remove spurious pole-zero pairs.

% Find negligible residues:
ii = find(abs(res) < cleanup_tol * norm(F, inf));
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
[pol, res, zer] = prz(r, z, f, w);

end % End of CLEANUP().


%% Automated choice of sample set

function [r, pol, res, zer, zj, fj, wj, errvec] = ...
    aaa_autoZ(F, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag)
%

% Flag if function has been resolved:
isResolved = 0;

% Main loop:
for n = 5:14
    % Sample points:
    % Next line enables us to do pretty well near poles
    Z = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), 1 + 2^n).';
    [r, pol, res, zer, zj, fj, wj, errvec] = aaa(F, Z, 'tol', tol, ...
        'mmax', mmax, 'cleanup', cleanup_flag, 'cleanuptol', cleanup_tol);
    
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
