function c = idlt(f)
%IDLT  Inverse discrete Legendre transform. 
%   V = DLT(C) returns a column vector V such that
%       V(k) = C(1)*P_0(x(k)) + C(2)*P_1(x(k)) + ... C(N)*P_{N-1}(x(k)), 
%   where P_j(x) is the degree j Legendre polynomial and x is the vector of
%   Gauss-Legendre nodes (as returned by x = legpts(size(C,1))).
%
% See also CHEBFUN.DLT, CHEBFUN.DCT.

% Nick Hale & Alex Townsend, Feb 2015.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: See N Hale & A Townsend - "A fast FFT-based DLT" for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(f, 1);

if ( N == 0 )
    
    % Trivial empty case.
    c = f;
    
elseif ( N < 5000 ) % <-- Determined expoerimentally.
    % Compute using direct (recurrence relation-based) method:
    c = idlt_direct(f);
    
else
    
    % Stage I:
    n = size(f, 1);
    [~, w] = legpts(n);                           % Legendre weights
    wf = bsxfun(@times, w.', f);                  % Scale by weights
    c = ndct_transpose(wf);                       % NDCT transpose
    
    % Stage II:
    c = leg2cheb(c, 'transpose');                 % LEG2CHEB transpose
    c = bsxfun(@times, (0:n-1).' + .5, c);        % Scaling
    
end

end

function v = idlt_direct(c)
%IDLT_DIRECT  Evaluate Legendre-Vandermonde matrix inverse times a vector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: Uses recurrence relation but does not form the matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(c, 1);
% Trivial case:
if ( N == 1 ), v = 1 + 0*c; return, end

% Compute the Legendre nodes:
[x, w] = legpts(N);
c = bsxfun(@times, w.', c);                   % Scale by weights
Pm1 = 1+0*x; P = x;                           % P_0 and P_1.
v = zeros(size(c));
v(1,:) = sum(c,1);
v(2,:) = x.'*c;
for n = 1:(N-2)                               % Recurrence relation:
    Pp1 = (2-1/(n+1))*(P.*x) - (1-1/(n+1))*Pm1;
    Pm1 = P; P = Pp1;
    v(n+2,:) = P.'*c;
end
v = bsxfun(@times, (0:N-1).' + .5, v);        % Scaling

end

function c = ndct_transpose(f)
%NDCT_TRANSPOSE  Evalute transpose of the NDCT operator.
%   See NDCT in DLT.M.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: Uses the 'cheb_1' version of the algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = eps;                                    % Tolerance
L = 18;                                       % Max number of terms in series

[n, m] = size(f); nn = (0:n-1).';             % Transform size
[~, ~, ~, t_leg] = legpts(n);                 % Legendre grid (in theta)
t_leg = t_leg(end:-1:1); f = f(end:-1:1,:);   % Get right-left ordering
t_cheb = ((0:n-1).'+.5)*pi./n;                % Chebyshev nodes (in theta)
dt = (t_leg - t_cheb);                        % (t_leg - t_cheb)
dt(n:-1:ceil((n+1)/2)) = -dt(1:ceil(n/2));    % More accurate near t = 0
dt = repmat(dt, 1, m); nn = repmat(nn, 1, m); % Support matrix input
NN = 1; scl = 1; sgn = -1;                    % Initialise
c = chebfun.dct( f, 2 );                      % ell = 0 term
for l = 1:(L-1)                               % Remaining terms in sum
    f = dt.*f; NN = nn.*NN; scl = scl*l;      % Update terms
    if ( mod(l, 2) )
        dc = (sgn/scl)*NN.*dst3_shifted_transpose(f); % Sine terms
    else
        dc = (sgn/scl)*NN.*chebfun.dct( f, 2 );       % Cosine terms
    end
    if ( l/2 == round(l/2) ),
        sgn = -sgn;                           % +--++--.
    end
    c = c + dc; 
    if ( norm(dc, inf) < tol ), break, end
end

% Ensure real/imag output for real/imag input:
if ( isreal(f) ), c = real(c); elseif ( isreal(1i*f) ), c = imag(c); end

end

function v = dst3_shifted_transpose(c)
% dst3_shifted_transpose is computed by a DST-III.
v = chebfun.dst( c, 2 );
v = [zeros(1, size(c, 2)) ; v(1:end-1,:)];
end
