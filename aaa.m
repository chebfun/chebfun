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
%   - 'cleanup', 'off': turns off automatic removal of numerical Froissart
%       doublets
%
%   One can also execute R = AAA(F), with no specification of a set Z.
%   This is equivalent to defining Z = LINSPACE(-1,1,1000) if F is a
%   function handle, Z = LINSPACE(A,B,1000) if F is a chebfun with
%   domain [A,B], and Z = LINSPACE(-1,1,LENGTH(F)) if F is a vector.
%
%   Reference:
%   [1] Yuji Nakatsukasa, Olivier Sete, Lloyd N. Trefethen, "The AAA algorithm
%   for rational approximation", arXiv:1612.00337.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% TODO:
% - Adaptive choice of the number of sample points when Z omitted and F is a
%   function handle or a chebfun.
% - Once adaptivity is implemented, add domains other than [-1, 1]:
%   possibility of other domains? 'unit' (roots of unity), [a,b].

% Parse inputs:
[F, Z, M, tol, mmax, cleanup_flag] = parseInputs(F, varargin{:});

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

% Construct function handle:
r = @(zz) reval(zz, zj, fj, wj);

% Compute poles, residues and zeros:
[pol, res, zer] = prz(r, zj, fj, wj);

if ( cleanup_flag )
    % Remove Froissart doublets:
    [r, pol, res, zer, zj, fj, wj] = cleanup(r, pol, res, zer, zj, fj, wj, Z, F);
end

end % of AAA()



%% parse Inputs:

function [F, Z, M, tol, mmax, cleanup_flag] = parseInputs(F, varargin)
% Input parsing for AAA.

% Sample points:
if ( isempty(varargin) || ~isnumeric(varargin{1}) )
    % Set of sample points Z is not given, default to equispaced points.
    numpts = 1000;
    if ( isnumeric(F) )
        % F is given as data values, pick same number of sample points:
        Z = linspace(-1, 1, length(F));
        
    elseif ( isa(F, 'chebfun') )
        if ( size(F, 2) ~= 1 )
            error('CHEBFUN:aaa:NonscalarF', 'The chebfun F must be scalar.')
        end
        % Default to the domain of F:
        dom = F.domain;
        Z = linspace(dom(1), dom(end), numpts);
        
    else
        % Default Z to [-1, 1]:
        Z = linspace(-1, 1, numpts);
    end
else %if ( isnumeric(varargin{1}) )
    Z = varargin{1};
    varargin(1) = [];
end

% Work with column vectors:
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
else
    error('CHEBFUN:aaa:UnknownF', 'Input for F not recognized.')
end

% Further parameters:

% Set defaults:
tol = 1e-13;        % Relative tolerance.
mmax = 100;         % Maximum number of terms.
cleanup_flag = 1;   % Cleanup on.

% Check further inputs:
while ( ~isempty(varargin) )
    if ( strncmpi(varargin{1}, 'tol', 3) )
        tol = varargin{2};
        varargin([1, 2]) = [];
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        mmax = varargin{2};
        varargin([1, 2]) = [];
    elseif ( strncmpi(varargin{1}, 'cleanup', 7) )
        if ( strncmpi(varargin{2}, 'off', 3) )
            cleanup_flag = 0;
        end
        varargin([1, 2]) = [];
    else
        error('CHEBFUN:aaa:UnknownArg', 'Argument unknown.')
    end
end

end

%% Evaluate rational function in barycentric form

function r = reval(zz, zj, fj, wj)
% Evaluate rational function in barycentric form.
zv = zz(:);                             % vectorize zz if necessary
CC = 1./bsxfun(@minus, zv, zj.');       % Cauchy matrix
r = (CC*(wj.*fj))./(CC*wj);             % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv)) = sum(wj.*fj)./sum(wj);

% Deal with NaN:
ii = find(isnan(r));
for jj = 1:length(ii)
    if ( isnan(zv(ii)) )
        % r(NaN) = NaN is fine.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj)) = fj(zv(ii(jj)) == zj);
    end
end

% Reshape to input format:
r = reshape(r, size(zz));

end % of REVAL()


%% Compute poles, residues and zeros.

function [pol, res, zer] = prz(r, zj, fj, wj)
% Compute poles, residues, and zeros of rational function in barycentric form.
m = length(wj);

% Compute poles via generalized eigenvalue problem:
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
% Remove zeros of denominator at infinity:
pol = pol(~isinf(pol));

% Compute residues via discretized Cauchy integral:
dz = 1e-5*exp(2i*pi*(1:4)/4);
res = r(bsxfun(@plus, pol, dz))*dz.'/4;

% Compute zeros via generalized eigenvalue problem:
E = [0 (wj.*fj).'; ones(m, 1) diag(zj)];
zer = eig(E, B);
% Remove zeros of numerator at infinity:
zer = zer(~isinf(zer));

end % end of PRZ()



%% Cleanup

function [r, pol, res, zer, z, f, w] = cleanup(r, pol, res, zer, z, f, w, Z, F)
% Remove spurious pole-zero pairs.

% Find negligible residues:
ii = find(abs(res) < 1e-13);
ni = length(ii);
if ni == 0
    % Nothing to do.
    return
elseif ni == 1
    fprintf('1 Froissart doublet.\n')
else
    fprintf('%d Froissart doublets.\n', ni)
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
    Z(Z == z(jj)) = [];
    F(Z == z(jj)) = [];
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

end % of CLEANUP()
