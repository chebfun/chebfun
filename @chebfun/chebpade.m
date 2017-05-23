function varargout = chebpade(F, m, n, varargin)
%CHEBPADE   Chebyshev-Pade approximation.
%   [P, Q, R_HANDLE] = CHEBPADE(F, M, N) computes polynomials P and Q of degree
%   M and N, respectively, such that the rational function P/Q is the type (M,
%   N) Chebyshev-Pade approximation of type Clenshaw-Lord to the CHEBFUN F. That
%   is, the Chebyshev series of P/Q coincides with that for the CHEBFUN F up to
%   the maximum possible order for the polynomial degrees permitted. R_HANDLE is
%   a function handle for evaluating the rational function.
%
%   [P, Q, R_HANDLE] = CHEBPADE(F, M, N, TYPE) allows one to additionally
%   specify the type of Chebyshev-Pade approximation sought. If TYPE is set to
%   'clenshawlord', the Clenshaw-Lord approximation as described above is used.
%   Alternatively, setting TYPE to 'maehly' computes a Maehly-type
%   approximation, which satisfies a linearized version of the Chebyshev-Pade
%   conditions.
%
%   [P, Q, R_HANDLE] = CHEBPADE(F, M, N, TYPE, K) uses only the K-th partial sum
%   in the Chebyshev expansion of F when computing the approximation. CHEPADE(F,
%   M, N, K) is shorthand for CHEBPADE(F, M, N, 'clenshawlord', K).
% 
%   In all of the above cases, if only one output argument is specified
%   then R_HANDLE is returned, while P and Q are returned if two output
%   arguments are specified. 
%
% See also AAA, CF, MINIMAX, PADEAPPROX, RATINTERP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Handle quasimatrices/array-valued CHEBFUNs.

% TODO:  References?

% Parse the inputs.
if ( nargin == 2 )     % CHEBPADE(F, M)
    M = -1;
    type = 'clenshawlord';
    n = 0;
elseif ( nargin == 3 ) % CHEBPADE(F, M, N)
    M = -1;
    type = 'clenshawlord';
    if ( ~isnumeric(n) )
        n = 0;
    end
elseif ( nargin == 4 ) % CHEBPADE(F, M, N, TYPE) or CHEPADE(F, M, N, K)
    if ( isnumeric(varargin{1}) )
        M = varargin{1};
        type = 'clenshawlord';
    elseif ( ischar(varargin{1}) )
        M = -1;
        type = varargin{1};
        if ( ~(strcmpi(type,'clenshawlord') || strcmpi(type,'maehly')) )
            error('CHEBFUN:CHEBFUN:chebpade:inputParameters', ...
                'Unrecognized sequence of input parameters.')
        end
    end
elseif ( nargin == 5 ) % CHEBPADE(F, M, N, TYPE, K) or CHEPADE(F, M, N, K, TYPE)
    if ( isnumeric(varargin{1}) && ischar(varargin{2}) )
        M = varargin{1};
        type = varargin{2};
    elseif ( ischar(varargin{1}) && isnumeric(varargin{2}) )
        M = varargin{2};
        type = varargin{1};
    else
        error('CHEBFUN:CHEBFUN:chebpade:inputParameters', ...
            'Unrecognized sequence of input parameters.')
    end
else
    error('CHEBFUN:CHEBFUN:chebpade:tooManyArgs', 'Too many arguments.')
end

if ( issing(F) )
    error('CHEBFUN:chebpade:singularFunction', ...
        'CHEBPADE does not currently support functions with singularities.');
end

% Compute the Chebyshev-Pade approximation of choice.
if ( strcmp(type,'clenshawlord') )
    [p, q, r_handle] = chebpadeClenshawLord(F, m, n, M);
elseif ( strcmp(type,'maehly') )
    [p, q, r_handle] = chebpadeMaehly(F, m, n);
else
    error('CHEBFUN:CHEBFUN:chebpade:type', ...
        'Approximation type must be either ''clenshawlord'' or ''maehly''.');
end

% Deal with transposed CHEBFUNs.
if ( F.isTransposed )
    p = p.';
    q = q.';
end

% Return the output:
outArgs = {p, q, r_handle};
if ( nargout <= 1 )
    varargout{1} = r_handle;
elseif ( nargout <= 3 )
    [varargout{1:nargout}] = outArgs{1:nargout};
else
    error('CHEBFUN:CHEBFUN:chebpade:nargout', ...
        'Incorrect number of output arguments.'); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clenshaw-Lord approximation.

function [p, q, r_handle] = chebpadeClenshawLord(F, m, n, M)

% Temporary degree variable in case m < n.
l = max(m, n);

% Limit to fixed number of coefficients if specified.
if ( M >= 0 )
    F = chebfun(@(x) feval(F, x), F.domain([1 end]), M + 1);
elseif ( numel(F.funs) > 1 )
    error('CHEBFUN:CHEBFUN:chebpade:chebpadeClenshawLord:multipleFuns', ...
        ['For a function with multiple funs, the number of coefficients ' ...
         'to be considered should be specified.']);
end

% Get the Chebyshev coefficients and pad if necessary.
c = chebcoeffs(F, length(F));
if ( length(F) < m + 2*n + 1 )
    %c = [c ; zeros(m + 2*n+1 - length(F), 1)];

    % Using random values appears to be more stable than using zeros?
    c = [c ; eps*randn(m + 2*n+1 - length(F), 1)];
    %rng('default')
end
c(1) = 2*c(1);

% Set up and solve Hankel system for denominator Laurent-Pade coefficients.
top = c(abs((m-n+1:m)) + 1); % Top row of Hankel system.
bot = c((m:m+n-1) + 1);      % Bottom row of Hankel system.
rhs = c((m+1:m+n) + 1);      % RHS of Hankel system.

beta = 1;
if ( n > 0 )
    beta = flipud([-hankel(top, bot)\rhs ; 1]);
end

% Use convolution to compute numerator Laurent-Pade coefficients.
c(1) = c(1)/2;
alpha = conv(c(1:l+1), beta);
alpha = alpha(1:l+1);
beta = beta.';

% Compute numerator Chebyshev-Pade coefficients.
pk = zeros(m + 1, 1);
D = zeros(l + 1, l + 1);
D(1:l+1,1:n+1) = alpha(:,ones(n + 1, 1)).*beta(ones(l + 1, 1),:);

pk(1) = sum(diag(D));
for k = 1:m
    pk(k+1) = sum([diag(D, k) ; diag(D, -k)]);
end

% Compute denominator Chebyshev-Pade coefficients.
qk = zeros(n + 1, 1);
for k = 1:n+1
    u = beta(1:n+2-k);
    v = beta(k:end);
    qk(k) = u*v.';
end

% Normalize the coefficients.
pk = pk/qk(1);
qk = 2*qk/qk(1);
qk(1) = 1;

% Form the outputs.
p = chebfun(pk, F.domain([1 end]), 'coeffs');
q = chebfun(qk, F.domain([1 end]), 'coeffs');
r_handle = @(x) feval(p, x)./feval(q, x);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maehly approximation.

function [p, q, r_handle] = chebpadeMaehly(F, m, n)

% Tolerance for determining negligible singular values.
tol = 1e-10;

% Get the Chebyshev coefficients and pad if necessary.
a = chebcoeffs(F, length(F));
if ( length(F) < m + 2*n + 1 )
    %a = [a ; zeros(m + 2*n + 1 - length(F), 1)];

    % Using random values appears to be more stable than using zeros?
    warning('CHEBFUN:CHEBFUN:chebpade:chebpadeMaehly:notEnoughCoeffs', ...
        ['Not enough coefficients given for type [' int2str(m) '/' ...
         int2str(n) '] approximation.  Assuming remainder are noise.']);
    a = [a ; eps*randn(m + 2*n + 1 - length(F), 1)];
    %rng('default')
end

% Form matrix for computing denominator coefficients.
row = (1:n);
col = (m+1:m+n).';
D = a(col(:,ones(n, 1)) + row(ones(n, 1),:) + 1) + ...
    a(abs(col(:,ones(n, 1)) - row(ones(n, 1),:)) + 1);
if ( n > m )
    D = D + a(1)*diag(ones(n - m, 1), m);
end

% Solve system for denominator coefficients.
if ( rank(D, tol) < min(size(D)) )
    % If system matrix is singular, reduce degrees first and try again.
    if ( m > 1 )
        [p, q, r_handle] = chebpade(F, m - 1, n, 'maehly');
        warning('CHEBFUN:CHEBFUN:chebpade:chebpadeMaehly:singularGoingLeft', ...
            ['Singular matrix encountered. Computing [' int2str(m - 1) ...
             '/' int2str(n) '] approximant.'])
    elseif ( n > 1 )
        [p, q, r_handle] = chebpade(F, m, n - 1, 'maehly');
        warning('CHEBFUN:CHEBFUN:chebpade:chebpadeMaehly:singularGoingUp', ...
            ['Singular matrix encountered. Computing [' int2str(m) ...
             '/' int2str(n - 1) '] approximant.'])
    else
        error('CHEBFUN:CHEBFUN:chebpade:chepadeMaehly:singularFail', ...
            ['Singular matrix encountered.  Cannot compute [1/1] ' ...
            'approximation.']);
    end
    return
else
    % Otherwise, solve for the denominator coefficients.
    qk = [1 ; -D\(2*a(m+2:m+n+1))];
end

% Compute numerator coefficients.
col = (1:m)';
B = a(col(:,ones(n, 1)) + row(ones(m, 1),:) + 1) + ...
    a(abs(col(:,ones(n, 1)) - row(ones(m, 1),:)) + 1);
mask = 1:(m + 1):min(m, n)*(m + 1);
B(mask) = B(mask) + a(1);

if ( m == 1 )
    B = B.';
end

B = [a(2:n+1).' ;  B];

if ( isempty(B) )
    pk = qk(1)*a(1:m+1);
else
    pk = .5*B*qk(2:n+1) + qk(1)*a(1:m+1);
end

% Form the outputs.
p = chebfun(pk, F.domain([1 end]), 'coeffs');
q = chebfun(qk, F.domain([1 end]), 'coeffs');
r_handle = @(x) feval(p, x)./feval(q, x);

end
