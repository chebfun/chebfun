function u = poisson( f, varargin )
%POISSON   Fast Poisson solver for the rectangle.
%   POISSON(F) solves laplacian(U) = F on the domain of F with zero
%   Dirichlet boundary conditions. That is, U satisfies
%
%     U_{x,x} + U_{y,y} = F, on [a,b]x[c,d], with U = 0 on boundary
%
%   The equation is solved using an adaptively determined discretization
%   size.
%
%   POISSON(F, G) solves using Dirichlet boundary conditions given by G. G
%   can be a scalar, a function handle, or any chebfun2 object satisfying
%   the Dirichlet data.
%
%   POISSON(F, G, N) is the same as POISSON(F, G), but uses an N x N tensor
%   product discretization to solve the equation.
%
%   POISSON(F, G, M, N) is the same as POISSON(F, G, N), but with an M x N
%   tensor product discretization, where N is the number of coeffcieints in
%   the x-direction and M is the number in the y-direction.
%
%   POISSON(F, G, N, METHOD) or POISSON(F, G, M, N, METHOD) is the same as
%   POISSON(F, G, N) or POISSON(F, G, M, N), respectively, except the
%   underlying matrix equation is solved with METHOD. Available methods
%   are:
%
%     'adi'            - alternating direction implicit method
%     'fadi'           - factored alternating direction implicit method
%     'bartelsStewart' - Bartels-Stewart algorithm
%
%   If METHOD is not supplied, then this command selects one (based on the
%   discretization size and the rank of the righthand side).
%
% EXAMPLE:
%   f = chebfun2( @(x,y) 1 + 0*x, [-1 2 0 1]);
%   u = chebfun2.poisson(f);
%   plot(u)
%
% See also DISKFUN/POISSON, SPHEREFUN/POISSON.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPER'S NOTE:
%
% METHOD: Spectral method (in coefficient space). We use a C^{(3/2)} basis
% to discretize the equation, resulting in a discretization of the form
% AX + XA = F, where A is a symmetric tridiagonal matrix.
%
% LINEAR ALGEBRA: Matrix equations. The matrix equation is solved by the
% alternating direction implicit (ADI) method.
%
% SOLVE COMPLEXITY:  O(M*N*log(MAX(M,N))*log(1/eps)) with M*N = total
% degrees of freedom.
%
% AUTHORS: Dan Fortunato (dan.fortunato@gmail.com)
%          Alex Townsend (townsend@cornell.edu)
%
% The fast Poisson solver is based on:
%
% D. Fortunato and A. Townsend, Fast Poisson solvers for spectral methods,
% Submitted, 2017.

% Solve for u on the same domain as f, adjust diffmat to include scaling:
dom = f.domain;
scl_x = (2/(dom(2)-dom(1)))^2;
scl_y = (2/(dom(4)-dom(3)))^2;

% Not enough input arguments so we call chebop2 to adaptively select a
% discretization size:
if ( nargin == 1 )

    % Call is POISSON(F) so set G = 0 and use chebop2 to determine the
    % discretization size:
    g = 0;
    N = chebop2.setupLaplace( f.domain );
    N.bc = g;
    u = N \ f;
    return

elseif ( nargin == 2 )

    % Call is POISSON(F, G) so use chebop2 to determine the discretization
    % size:
    g = varargin{1};
    if ( ~isa(g, 'chebfun2') )
        g = chebfun2(g, f.domain);
    end
    N = chebop2.setupLaplace( f.domain );
    % Note: chebfun2/subsref (e.g. g(-1,:)) does not work here, so we use
    %       feval instead. Is this a bug?
    N.lbc = feval(g, f.domain(1), ':');
    N.rbc = feval(g, f.domain(2), ':');
    N.dbc = feval(g, ':', f.domain(3));
    N.ubc = feval(g, ':', f.domain(4));
    u = N \ f;
    return

end

% We are given a discretization size. It's solve time!
g = varargin{1};
m = varargin{2};
method = '';
if ( nargin == 3 )
    % Call must be POISSON(F, G, N) so employ an NxN discretization:
    n = m; % square discretization
elseif ( nargin == 4 )

    % The user typed one of the following:
    % POISSON(F, G, N, METHOD) or POISSON(F, G, M, N).

    if ( ischar( varargin{3} ) )
        % Call is POISSON(F, G, N, METHOD):
        method = varargin{3};
        n = m; % square discretization
    else
        % Call is POISSON(F, G, M, N):
        n = varargin{3};
    end

elseif ( nargin == 5 )

    % We are given a discretization size and a method. It's solve time!
    % Call must be POISSON(F, G, M, N, METHOD).
    n = varargin{3};
    method = varargin{4};

else
    error('CHEBFUN2:POISSON:NARGIN', ...
        'Too many input arguments to chebfun2.poisson().');
end

% Set the error tolerance for solve:
tol = chebfun2eps();

% Compute the Chebyshev coefficients of rhs:
[Cf, Df, Rf] = coeffs2(f, m, n);

% Solver only deals with zero homogeneous Dirichlet conditions. Therefore,
% if nonzero Dirichlet conditions are given, we solve lap(u) = f with u|bc = g
% as u = v + w, where v|bc = g, and lap(w) = f - lap(v), w|bc = 0:
if ( isa(g, 'double') || isa(g, 'chebfun2') || isa(g, 'function_handle') )

    % Make double or function handle into chebfun2:
    if ( isa(g, 'double') || isa(g, 'function_handle') )
        g = chebfun2(g, f.domain);
    end

    if ( g.domain == f.domain )
        % Adjust the rhs, if nonzero:
        lapg = lap(g);
        if ( ~iszero( lapg ) )
            [Cg, Dg, Rg] = coeffs2(lapg, m, n);
            Cf = [Cf Cg];
            Z = zeros(size(Df,1),size(Dg,1));
            Df = [Df Z ; Z' -Dg];
            Rf = [Rf Rg];
        end
    else
        error('CHEBFUN2:POISSON:BC', ...
            'Dirichlet data should be on the same domain as F.');
    end

else
    error('CHEBFUN2:POISSON', ...
        'Dirichlet data needs to be given as a scalar or function.')
end

% Convert rhs to C^{(3/2)} coefficients:
Cf = cheb2ultra( Cf );
Rf = cheb2ultra( Rf );

% Construct M, the multiplication matrix for (1-x^2) in the C^(3/2) basis
jj = (0:n-1)';
dsub = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj+2);
dsup = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj);
d = -dsub - dsup;
Mn = spdiags([dsub d dsup], [-2 0 2], n, n);
% Construct D^{-1}, which undoes the scaling from the Laplacian identity
invDn = spdiags(-1./(jj.*(jj+3)+2), 0, n, n);
Tn = scl_y * invDn * Mn;

jj = (0:m-1)';
dsub = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj+2);
dsup = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj);
d = -dsub - dsup;
Mm = spdiags([dsub d dsup], [-2 0 2], m, m);
invDm = spdiags(-1./(jj.*(jj+3)+2), 0, m, m);

% Construct T = D^{-1} * M:
Tm = scl_x * invDm * Mm;
Cf = invDm * Cf;
Rf = invDn * Rf;

switch lower(method)

    case 'bartelsstewart'

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%  Bartels-Stewart method %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve TmX + XTn' = F using Bartels-Stewart, which requires O(n^3)
        % operations:

        X = chebop2.bartelsStewart(Tm, eye(n), eye(m), Tn, Cf*Df*Rf.', 0, 0);

        % Convert back to Chebyshev
        X = ultra1mx2cheb( ultra1mx2cheb( X ).' ).';
        u = chebfun2( X, f.domain, 'coeffs' );

    case {'adi', 'fadi', ''}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%  Alternating Direction Implicit method %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve TmX + XTn' = F using ADI, which requires
        % O(n^2log(n)log(1/eps)) operations:

        % An ADI method will be used (either given by the user, or selected
        % by us.)

        % Compute ADI shifts
        a =  -4/pi^2 * scl_y;
        b = -39*n^-4 * scl_y;
        c =  39*m^-4 * scl_x;
        d =   4/pi^2 * scl_x;
        [p, q] = chebop2.adiShifts(a, b, c, d, tol);

        if ( isempty(method) )
            % Let's go and pick a good method to use:
            % Test if we should use ADI or FADI:
            rho = size(Cf,2); % Rank of rhs
            adi_test = ( min(m,n) < rho*numel(p)/2 ); % Worth doing FADI?
            if ( adi_test )
                method = 'adi';
            else
                method = 'fadi';
            end
        end

        % Solve matrix equation:
        if ( strcmpi(method, 'adi') )
            % Run the ADI method:
            X = chebop2.adi(Tm, -Tn', Cf*Df*Rf.', p, q );

            % Convert back to Chebyshev
            X = ultra1mx2cheb( ultra1mx2cheb( X ).' ).';
            u = chebfun2( X, f.domain, 'coeffs' );

        else
            % Run the FADI method:
            [UX, DX, VX] = chebop2.fadi(Tm, -Tn, Cf*Df, Rf, p, q);

            % Convert back to Chebyshev:
            UX = ultra1mx2cheb(UX);
            VX = ultra1mx2cheb(VX);

            UX = chebfun(UX, dom(3:4), 'coeffs');
            VX = chebfun(VX, dom(1:2), 'coeffs');
            u = chebfun2();
            u.cols = UX;
            u.rows = VX;
            u.pivotValues = 1./diag(DX);
            u.domain = dom;
        end

    otherwise
        error('CHEBFUN2:POISSON:SOLVER', ...
            'Method supplied to chebfun2.poisson() is not recognized.');
end

% Add back in the boundary data:
u = u + g;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CONVERSION CODES %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = cheb2ultra( X )
% CHEB2ULTRA   Convert vector of Chebyshev coefficients to C^(3/2).
%
%    CHEB2ULTRA(X) applies the conversion to each column of X if X is a
%    matrix.

% First convert the matrix of Chebyshev coefficients to a matrix of
% Legendre coefficients:
m = size( X, 1 );
if ( m <= 10000 ) % Determined experimentally
    S = cheb2leg_mat( m );
    X = S * X;
else
    X = cheb2leg( X );
end

% Now, convert the matrix of Legendre coefficients to a matrix of
% ultraspherical coefficients:
S = leg2ultra_mat( m );
X = S * X;

end

function X = ultra1mx2cheb( X )
% ULTRA1MX2CHEB    Convert vector of (1-x^2)C^(3/2) coefficients to
% Chebyshev.
%
%  ULTRA1MX2CHEB(X) applies the conversion each column of X if X is
%    matrix.

% First, convert the matrix of (1-x^2)C^(3/2)(x) coefficients
% to Legendre coefficients:
m = size( X, 1 );

S = ultra1mx2leg_mat( m );
X = S * X;

% Now, convert the matrix of Legendre coefficient to a matrix of Chebyshev
% coefficients:
if ( m <= 10000 ) % Determined experimentally
    S = leg2cheb_mat( m );
    X = S * X;
else
    X = leg2cheb( X );
end

end

function S = leg2ultra_mat( n )
% LEG2ULTRA_MAT Conversion matrix from Legendre coefficients to C^(3/2).
%
% Given coefficients in the Legendre basis the C^(3/2) coefficients
% can be computed via
%
%     c = rand(10, 1);    % Legendre coefficients
%     S = leg2ultra_mat( length(c) ); % conversion matrix
%     d = S * c;           % C^(3/2) coefficients

% Alex Townsend, 5th May 2016

lam = 1/2;
dg = lam./(lam + (2:n-1))';
v  = [1 ; lam./(lam+1) ; dg];
w  = [0 ; 0 ; -dg];
S  = spdiags( [v w], [0 2], n, n );

end

function S = ultra1mx2leg_mat( n )
% ULTRA1MX2LEG_MAT Conversion matrix for (1-x^2)C^(3/2) to Legendre.
%
% Given coefficients in the (1-x^2)C^(3/2) basis the Legendre coefficients
% can be computed via
%
%     c = rand(10, 1);     % (1-x^2)C^(3/2) coefficients
%     S = ultra1mx2leg_mat( length(c) ); % conversion matrix
%     c_leg = S * c;       % Legendre coefficients
%

% Alex Townsend, 5th May 2016

d = ones(n, 1);
S = spdiags(((1:n).*(2:(n+1))./2./(3/2:n+1/2))', 0, n, n);
S = spdiags( [d,-d], [0,-2], n, n ) * S;

end

function L = cheb2leg_mat( N )
% CHEB2LEG_MAT Construct the cheb2leg conversion matrix.

% This for-loop is a faster and more accurate way of doing:
% Lambda = @(z) exp(gammaln(z+1/2) - gammaln(z+1));
% vals = Lambda( (0:2*N-1)'/2 );
vals = zeros(2*N,1);
vals(1) = sqrt(pi);
vals(2) = 2/vals(1);
for i = 2:2:2*(N-1)
    vals(i+1) = vals(i-1)*(1-1/i);
    vals(i+2) = vals(i)*(1-1/(i+1));
end

L = zeros(N, N);
for j = 0:N-1
    for k = j+2:2:N-1
        L(j+1, k+1) = -k*(j+.5)*(vals((k-j-2)+1)./(k-j)).*(vals((k+j-1)+1)./(j+k+1));
    end
end
c = sqrt(pi)/2;
for j = 1:N-1
    L(j+1, j+1) = c./vals( 2*j+1 );
end
L(1,1) = 1;

end

function M = leg2cheb_mat( N )
% LEG2CHEB_MAT Construct the leg2cheb conversion matrix.

% This for-loop is a faster and more accurate way of doing:
% Lambda = @(z) exp(gammaln(z+1/2) - gammaln(z+1));
% vals = Lambda( (0:2*N-1)'/2 );
vals = zeros(2*N,1);
vals(1) = sqrt(pi);
vals(2) = 2/vals(1);
for i = 2:2:2*(N-1)
    vals(i+1) = vals(i-1)*(1-1/i);
    vals(i+2) = vals(i)*(1-1/(i+1));
end

M = zeros(N, N);
for j = 0:N-1
    for k = j:2:N-1
        M(j+1, k+1) = 2/pi*vals((k-j)+1).*vals((k+j)+1);
    end
end
M(1,:) = .5*M(1,:);

end