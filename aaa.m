function [r, pol, res, zer, zj, fj, wj, errvec, wt] = aaa(F, varargin)
%AAA   AAA and AAA-Lawson (near-minimax) real or complex rational approximation.
%   R = AAA(F, Z) computes the AAA rational approximant R (function handle) to
%   data F on the set of sample points Z.  F may be given by its values at Z,
%   or as a function handle or chebfun.  R = AAA(F, Z, 'degree', N) attempts to
%   compute the minimax approximation of degree N (i.e., rational type (N,N)).
%
%   [R, POL, RES, ZER] = AAA(F, Z) returns vectors of poles, residues, and zeros
%   of R.
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ] = AAA(F, Z) also returns the vectors of
%   support points ZJ, approximation values FJ = r(ZJ), and weights WJ of the
%   barycentric representation of R. 
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC] = AAA(F, Z) also returns the vector
%   of errors ||f-r||_infty in successive iterative steps of AAA.  Note that the
%   rational degrees are not 1,...,length(ERRVEC) but 0,...,length(ERRVEC)-1.
%
%   R = AAA(F, Z, NAME, VALUE) sets the following parameters:
%   - 'tol', TOL: relative tolerance (default TOL = 1e-13),
%   - 'degree', N: maximal degree (default N = 99). 
%       The output rational approximant will be of degree at most N.  Like
%       'mmax', N+1, except that Lawson is turned on by default: see below.
%   - 'mmax', MMAX: maximal number of terms in the barycentric representation
%       (default MMAX = 100). R will be of degree at most MMAX-1.  Like
%       'degree', MMAX-1, except that Lawson is not turned on by default.
%   - 'dom', DOM: domain (default DOM = [-1, 1]).  No effect if Z is provided.
%   - 'cleanup', 'off' or 0: turns off automatic removal of Froissart doublets.
%   - 'cleanuptol', CLEANUPTOL: cleanup tolerance (default CLEANUPTOL = TOL).
%       Poles with residues less than this number times the geometric mean size
%       of F times the minimum distance to Z are deemed spurious by the cleanup
%       procedure.  If TOL = 0, then CLEANUPTOL defaults to 1e-13.
%   - 'lawson', NLAWSON: take NLAWSON iteratively reweighted least-squares steps
%       to bring approximation closer to minimax.  Specifying NLAWSON = 0 
%       ensures there is no Lawson iteration.  See next paragraph.
%
%   If 'degree' is specified and 'lawson' is not, AAA attempts to find a minimax
%   approximant of degree N by AAA-Lawson iteration.  This will generally be
%   successful only if the minimax error is well above machine precision, and
%   is more reliable for complex problems than real ones.  If 'degree' and 
%   'lawson' are both specified, then exactly NLAWSON Lawson steps are taken
%   (so NLAWSON = 0 corresponds to AAA approximation with no Lawson iteration).
%   The final weight vector WT of the Lawson iteration is available with
%   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC, WT] = AAA(F, Z).
%
%   Note that R may have fewer than N poles and zeros.  This may happen, for
%   example, if N is too large, or if F is even and N is odd, or if F is odd
%   and N is even.
%
%   One can also execute R = AAA(F), with no specification of a set Z.  If F is
%   a vector, this is equivalent to R = AAA(F, Z) with
%   Z = LINSPACE(-1, 1, LENGTH(F)).  If F is a function handle, AAA attempts
%   to resolve F on [-1,1] by default.
%
%   This standalone code works in GNU Octave as well as MATLAB.
%
% Examples:
%   r = aaa(@exp); xx = linspace(-1,1); plot(xx,r(xx)-exp(xx))
%
%   r = aaa(@exp,'degree',4); xx = linspace(-1,1); plot(xx,r(xx)-exp(xx))
%
%   X = linspace(-1,1,30); r = aaa(gamma(X),X);
%   fplot(r,[-5,5]), axis([-5 5 -15 15]), grid on 
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
%   References on AAA and AAA-Lawson, respectively:
%
%   [1] Y. Nakatsukasa, O. Sete, and L. N. Trefethen, "The AAA algorithm
%   for rational approximation", SIAM J. Sci. Comp. 40 (2018), A1494-A1522.
%
%   [2] Y. Nakatsukasa and L. N. Trefethen, An algorithm for real and
%   complex rational minimax approximation, SIAM J. Sci. Comp. 42 (2020),
%   A3157-A3179.
%
% See also AAATRIG, CF, CHEBPADE, MINIMAX, PADEAPPROX, RATINTERP.

% Copyright 2023 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
[F, Z, M, dom, tol, mmax, cleanup_flag, cleanup_tol, needZ, mmax_flag, ...
    nlawson, degree_flag, degree] = parseInputs(F, varargin{:});

if ( needZ )
    % Z was not provided.  Try to resolve F on its domain.
    [r, pol, res, zer, zj, fj, wj, errvec] = ...
        aaa_autoZ(F, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag, ...
            nlawson, degree_flag, degree);
    return
end

% Remove infinite or NaN function values and repeated entries:
toKeep = ~isinf(F) & ~isnan(F);
F = F(toKeep); Z = Z(toKeep);
[Z, uni] = unique(Z,'stable'); F = F(uni);

% Initialization for AAA iteration:
M = length(Z);
abstol = tol*norm(F, inf);                 % Absolute tolerance
J = (1:M)';
zj = []; fj = []; C = []; A = [];
errvec = [];
R = mean(F)*ones(size(J));

% AAA iteration:
for m = 1:mmax
    % Introduce next support point:
    [~, jj] = max(abs(F(J) - R(J)));       % Select next support point
    zj = [zj; Z(J(jj))];                   % Update support points
    fj = [fj; F(J(jj))];                   % Update data values
    C = [C 1./(Z - Z(J(jj)))];             % Next column of Cauchy matrix
    J(jj) = [];                            % Update index vector
    A = [A, (F-fj(end)).*C(:,end)];        % Update Loewner matrix

    % Compute weights:
    if ( length(J) >= m )                  % The usual tall-skinny case
        [~, S, V] = svd(A(J,:), 0);        % Reduced SVD
        s = diag(S);
        mm = find( s == min(s) );          % Treat case of multiple min sing val
        nm = length(mm);
        wj = V(:,mm)*ones(nm,1)/sqrt(nm);  % Aim for non-sparse wt vector
    elseif ( length(J) >= 1 )
        V = null(A(J,:));                  % Fewer rows than columns
        nm = size(V,2);                    
        wj = V*ones(nm,1)/sqrt(nm);        % Aim for non-sparse wt vector
    else
        wj = ones(m,1)/sqrt(m);            % No rows at all (needed for Octave)
    end
    
    % Compute rational approximant:
    i0 = find(wj~=0);                      % Omit columns with wj = 0
    N = C(:,i0)*(wj(i0).*fj(i0));          % Numerator
    D = C(:,i0)*wj(i0);                    % Denominator
    R = N./D;
    Dinf = isinf(D);
    R(Dinf) = F(Dinf);                     % Interpolate at supp pts with wj~=0
    
    % Check if converged:
    maxerr = norm(F - R, inf);
    errvec = [errvec; maxerr];
    if ( maxerr <= abstol )
        break
    end
end
maxerrAAA = maxerr;                        % Error at end of AAA 

% We now enter Lawson iteration: barycentric IRLS = iteratively reweighted
% least-squares if 'lawson' is specified with NLAWSON > 0 or 'mmax' is
% specified and 'lawson' is not.  In the latter case the number of steps
% is chosen adaptively.  Note that the Lawson iteration is unlikely to be
% successful when the errors are close to machine precision.

wj0 = wj; fj0 = fj;                        % Save params in case Lawson fails
wt = NaN(M,1); wt_new = ones(M,1);
if ( nlawson > 0 )                         % Lawson iteration

    maxerrold = maxerrAAA;
    maxerr = maxerrold;
    nj = length(zj);
    A = [];
    for j = 1:nj                           % Cauchy/Loewner matrix
        A = [A 1./(Z-zj(j)) F./(Z-zj(j))];
    end
    for j = 1:nj
        [i,~] = find(Z==zj(j));            % support pt rows are special
        A(i,:) = 0;
        A(i,2*j-1) = 1;
        A(i,2*j) = F(i);
    end
    stepno = 0;
    while ( (nlawson < inf) && (stepno < nlawson) ) || ...
          ( (nlawson == inf) && (stepno < 20) ) || ...
          ( (nlawson == inf) && (maxerr/maxerrold < .999) && (stepno < 1000) ) 
        stepno = stepno + 1;
        wt = wt_new;
        W = spdiags(sqrt(wt),0,M,M);
        [~,~,V] = svd(W*A,0);
        c = V(:,end);
        denom = zeros(M,1); num = zeros(M,1);
        for j = 1:nj
            denom = denom + c(2*j)./(Z-zj(j));
            num = num - c(2*j-1)./(Z-zj(j));
        end
        R = num./denom;
        for j = 1:nj
            [i,~] = find(Z==zj(j));        % support pt rows are special
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
    if ( (maxerr > maxerrAAA) && (nlawson == Inf) )
        wj = wj0; fj = fj0; 
    end
end

% Remove support points with zero weight:
I = find(wj == 0);
zj(I) = []; wj(I) = []; fj(I) = [];

% Construct function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, zj, fj, wj);
[pol, res, zer] = prz(zj, fj, wj);

if ( cleanup_flag == 1 && nlawson == 0 )      % Remove Froissart doublets
    [r, pol, res, zer, zj, fj, wj] = ...
        cleanup(r, pol, res, zer, zj, fj, wj, Z, F, cleanup_tol);
elseif ( cleanup_flag == 2 && nlawson == 0 )  % Alternative cleanup.  Currently
    a.zj = zj; a.fj = fj; a.wj = wj;          % an undocumented feature,
    a.Z = Z; a.F = F;                         % pending further investigation.
    a.cleanup_tol = max(cleanup_tol, eps);
    c = cleanup2(a);
    zj = c.zj; fj = c.fj; wj = c.wj;
    r = @(zz) reval(zz, zj, fj, wj);
    [pol, res, zer] = prz(zj, fj, wj);
end

end % of AAA()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARSEINPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, Z, M, dom, tol, mmax, cleanup_flag, cleanup_tol, ...
    needZ, mmax_flag, nlawson, degree_flag, degree] = parseInputs(F, varargin)

% Check if F is empty:
if ( isempty(F) )
    error('AAA:emptyF', 'No function given.')
elseif ( isa(F, 'chebfun') )
    if ( size(F, 2) ~= 1 )
        error('AAA:nColF', 'Input chebfun must have one column.')
    end
end

% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    % Z is given.
    Z = varargin{1};
    if ( isempty(Z) )
        error('AAA:emptyZ', ...
            'If sample set is provided, it must be nonempty.')
    end
    varargin(1) = [];
end

% Set defaults for other parameters:
tol = 1e-13;                   % Relative tolerance
mmax = 100;                    % Maximum number of terms
degree = NaN;                  % Specified degree
cleanup_tol = 1e-13;           % Cleanup tolerance
nlawson = Inf;                 % Number of Lawson steps (Inf means adaptive)
% Domain:
if ( isa(F, 'chebfun') )
    dom = F.domain([1, end]);
else
    dom = [-1, 1];
end
cleanup_flag = 1;              % Cleanup on
mmax_flag = 0;                 % Checks if mmax manually specified
degree_flag = 0;               % Checks if degree specified
cleanup_set = 0;               % Checks if cleanup_tol manually specified
while ( ~isempty(varargin) )   % Check if parameters have been provided
    if ( strncmpi(varargin{1}, 'tol', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tol = varargin{2};
            if ( ~cleanup_set && tol > 0 ) % If not set, set cleanup_tol to tol
              cleanup_tol = tol;
            end
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'degree', 6) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            if ( mmax_flag == 1 ) && ( mmax ~= varargin{2}+1 )
                error('AAA:degmmaxmismatch', ' mmax must equal degree+1.')
            end            
            degree = varargin{2};
            mmax = degree + 1;
            mmax_flag = 1; 
            degree_flag = 1;
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )            
            if ( mmax_flag == 1 ) && ( mmax ~= varargin{2})                
                error('AAA:degmmaxmismatch', ' mmax must equal degree+1.')
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
                warning('AAA:dom', ...
                    ['Given domain does not match that of the chebfun.\n', ...
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
        error('AAA:UnknownArg', 'Argument unknown.')
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
            error('AAA:lengthFZ', ...
                'Inputs F and Z must have the same length.')
        end
    elseif ( ischar(F) )
        % F is given as a string input. Convert it to a function handle.
        F = inline(vectorize(F));
        F = F(Z);
    else
        error('AAA:UnknownF', 'Input for F not recognized.')
    end
    
else
    % Z was not given.  Set flag that Z needs to be determined.
    % Also set Z and M since they are needed as output.
    needZ = 1;
    Z = [];
    M = length(Z);
end

if ( ~degree_flag && (nlawson == Inf) )
    nlawson = 0;               
end

end % End of PARSEINPUTS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CLEANUP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
for j = 1:length(Zdistances)
    Zdistances(j) = min(abs(pol(j)-Z));
end
ii = find(abs(res)./Zdistances < cleanup_tol * geometric_mean_of_absF);
ni = length(ii);
if ( ni == 0 )
    return
elseif ( ni == 1 )
    warning('AAA:Froissart','1 Froissart doublet');
else
    warning('AAA:Froissart',[int2str(ni) ' Froissart doublets']);
end

% For each spurious pole find and remove closest support point:
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
    
    % Remove support point(s):
    z(jj) = []; f(jj) = [];
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
C = 1./(Z-z.');               % Cauchy matrix.
A = SF*C - C*Sf;              % Loewner matrix.

% Solve least-squares problem to obtain weights:
[~, ~, V] = svd(A, 0);
w = V(:,m);

% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, z, f, w);
[pol, res, zer] = prz(z, f, w);

end % End of CLEANUP.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CLEANUP2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        warning('AAA:Froissart',...
            ['1 Froissart doublet, niter = ', int2str(niter)]);
    else
        warning('AAA:Froissart',...
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
    C = 1./(Z-z.');             % Cauchy matrix.
    A = SF*C - C*Sf;            % Loewner matrix.
    
    % Solve least-squares problem to obtain weights:
    [~, ~, V] = svd(A, 0);
    w = V(:,m);
    
    % Compute poles, residues and zeros:
    [pol, res, zer] = prz(z, f, w);
end % End of while loop

c.zj = z; c.fj = f; c.wj = w;

end  % End of CLEANUP2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   AAA_AUTOZ   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, pol, res, zer, zj, fj, wj, errvec] = ...
    aaa_autoZ(F, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag, ...
               nlawson, degree_flag, degree)
% Automated choice of sample set

% Flag if function has been resolved:
isResolved = 0;

% Main loop:
for n = 5:14
    % Sample points:
    % Next line enables us to do pretty well near poles
    Z = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), 1 + 2^n).';
    if degree_flag
       [r, pol, res, zer, zj, fj, wj, errvec] = aaa(F, Z, 'tol', tol, ...
          'mmax', mmax, 'cleanup', cleanup_flag, 'cleanuptol', cleanup_tol, ...
          'lawson', nlawson, 'degree', degree);
    else
       [r, pol, res, zer, zj, fj, wj, errvec] = aaa(F, Z, 'tol', tol, ...
          'mmax', mmax, 'cleanup', cleanup_flag, 'cleanuptol', cleanup_tol, ...
          'lawson', nlawson);
    end
    % Test if rational approximant is accurate:
    abstol = tol * norm(F(Z), inf);
    
    % On Z(n):
    err(1,1) = norm(F(Z) - r(Z), inf);
    
    Zrefined = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), ...
        round(1.5 * (1 + 2^(n+1)))).';
    err(2,1) = norm(F(Zrefined) - r(Zrefined), inf);
    if ( all(err < abstol) )
        % Final check that the function is resolved, inspired by sampleTest().
        % Pseudo random sample points in [-1, 1]:
        xeval = [-0.357998918959666; 0.036785641195074];
        % Scale to dom:
        xeval = (dom(2) - dom(1))/2 * xeval + (dom(2) + dom(1))/2;
        
        if ( norm(F(xeval) - r(xeval), inf) < abstol )
            isResolved = 1;
            break
        end
    end
end

if ( ( isResolved == 0 ) && ~mmax_flag )
    warning('AAA:notResolved', ...
        'Function not resolved using %d pts.', length(Z))
end

end % End of AAA_AUTOZ.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PRZ   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pol, res, zer] = prz(zj, fj, wj)
%   Compute poles, residues, and zeros of rational fun in barycentric form.

% Compute poles via generalized eigenvalue problem:
m = length(wj);
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
pol = pol(~isinf(pol));

% Compute residues via formula for res of quotient of analytic functions:
N = @(t) (1./(t-zj.')) * (fj.*wj);
Ddiff = @(t) -((1./(t-zj.')).^2) * wj;
res = N(pol)./Ddiff(pol);

% Compute zeros via generalized eigenvalue problem:
E = [0 (wj.*fj).'; ones(m, 1) diag(zj)];
zer = eig(E, B);
zer = zer(~isinf(zer));

end % End of PRZ.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   REVAL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = reval(zz, zj, fj, wj)
%   Construct function handle to evaluate rational function in barycentric form.

zv = zz(:);                         % vectorize zz if necessary
CC = 1./(zv-zj.');                  % Cauchy matrix
r = (CC*(wj.*fj))./(CC*wj);         % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv)) = sum(wj.*fj)./sum(wj);

% Deal with NaN:
ii = find(isnan(r));
for jj = 1:length(ii)
    if ( isnan(zv(ii(jj))) || ~any(zv(ii(jj)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj)) = fj(zv(ii(jj)) == zj);
    end
end

% Reshape to input format:
r = reshape(r, size(zz));

end % End of REVAL.
