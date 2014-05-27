function c_leg = cheb2leg(c_cheb, M)
%LEG2CHEB convert Legendre coefficients to Chebyshev coefficients. 
%   C_LEG = CHEB2LEG(C_CHEB) converts the vector C_CHEB of Chebyshev
%   coefficients to a vector C_CHEB of Legendre coefficients such that
%   C_CHEB(N)*T0 + ... + C_CHEB(1)*T{N-1} = C_LEG(N)*P0 + ... + C_LEG(1)*P{N-1}.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This algorithm requires O( N(log N)^2 / log log N) operations and is based on
% rewritting an asymptotic formula for Legendre polynomials in a way that can be
% evaluated using discrete cosine transforms. For more details see:
%   N. Hale and A. Townsend, A fast, simple, and stable Chebyshev-Legendre
%   transform using an asymptotic formula, SISC (accepted) 2013.

c_cheb = c_cheb(:);                           % Make column vector.
c_cheb = flipud(c_cheb);                      % Lowest order coeffs first.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin == 1), M = 10; end                 % No. of terms in expansion.
N = length(c_cheb) - 1; NN = (0:N).';         % Degree of polynomial.
nM0 = min(floor(.5*(.25*eps*pi^1.5*gamma(M+1)/gamma(M+.5)^2)^(-1/(M+.5))), N);
aM = min(1/log(N/nM0), .5);                   % Block reduction factor (alpha_M)
K = ceil(log(N/nM0)/log(1/aM));               % Number of block partitions

% Use direct approach if N is small:
if ( M == 0 || N < 513 || K == 0 ), c_leg = cheb2leg_direct(c_cheb); return, end

f = dct1([c_cheb ; zeros(N,1)]);              % Values on a 2*N+1 Cheb grid.
wf = cheb2_quadwts(2*N+1).'.*f;               % Scale f by C-C weights.
t = pi*(0:2*N)'/(2*N);                        % 2*N+1 theta grid.
nM = ceil(aM.^(K-1:-1:0)*N);                  % n_M for each block.
jK = zeros(K, 2);                             % Block locations in theta.
for k = 1:K % Find where curve intersects a^k*N:
    tmp = find(t >= asin(nM0./nM(k)), 1) - 4; % Where curve intersects aM^k*N.    
    jK(k,:) = [tmp+1, 2*N+1-tmp];             % Collect indicies.
end
C = constantOutTheFront(N);                   % Scaling in asymptotic expansion.
nM(end) = N+2; % For convenience (avoids treating final block differently).

%% %%%%%%%%%%%%%%%%%%% Recurrence / boundary region %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_rec = zeros(N+1, 1);
jK2 = [jK(1,2)+1, 2*N+1 ; jK(1:K-1,:) ; 1, 2*N+1];
for k = 1:K % Loop over the block partitions:
    j_bdy = [jK2(k+1,1):jK2(k,1)-1, jK2(k,2)+1:jK2(k+1,2)]; % Boundary indicies.
    x_bdy = cos(t(j_bdy));                             % Boundary x values.
    wf_bdy = wf(j_bdy);                                % w.*f at boundary.
    tmp = zeros(N+1,1);                                % Local LHS storage.
    tmp(1:2) = [sum(wf_bdy) ; x_bdy'*wf_bdy];          % First two terms.
    Pm2 = 1; Pm1 = x_bdy;                              % Initialise recurrence.
    for n = 1:min(nM(k)-3, N-1) % Recurrence:
        P = (2-1/(n+1))*Pm1.*x_bdy - (1-1/(n+1))*Pm2; 
        Pm2 = Pm1; Pm1 = P;
        tmp(n+2) = P'*wf_bdy;                          % Update local LHS.
    end
    c_rec = c_rec + tmp;                               % Global correction LHS.
end

%% %%%%%%%%%%%%%%%%% Asymptotics / interior region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_leg = zeros(N+1, 1);                                 % Initialise output.
dst1Transpose([], 1);                                  % Clear persistent storage.
for k = 1:K-1 % Loop over the block partitions.
    c_k = zeros(N+1, 1);                               % Initialise local LHS.
    hm = ones(N+1,1); hm([1:nM(k)-1,nM(k+1):end]) = 0; % Initialise h_m.
    j_k = jK(k,1):jK(k,2);                             % t indices of kth block.
    t_k = pi/2 + 0*t; t_k(j_k) = t(j_k);               % t indices of kth block.
    wf_k = 0*wf; wf_k(j_k) = wf(j_k);                  % RHS w.*f of kth block.
    denom = 1./sqrt(2*sin(t_k));                       % Initialise denominator.
    for m = 0:M-1 % Terms in the asymptotic formula:
        denom = (2*sin(t_k)).*denom;                   % Update denominator.
        u = sin((m+.5)*(.5*pi-t_k))./denom;            % Trig terms:
        v = cos((m+.5)*(.5*pi-t_k))./denom;
        Tv = dct1(wf_k.*v);                            % Compute T'v.
        Uu = dst1Transpose(wf_k.*u);                   % Compute U'u*sin(t).
        c_k = c_k + hm.*([0 ; Uu(1:N)] + Tv(1:N+1));   % Update LHS vector.
        hm = ((m+0.5)^2./((m+1)*(NN+m+1.5))).*hm;      % Update h_m.
    end
    c_leg = c_leg + C.*c_k;                            % Append to global LHS.
end
dst1Transpose([], 1);                                  % Clear persistent storage.

%% %%%%%%%%%%%%%%%%%%%%%%%%% Combine for result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = (2*(0:N).'+1)/2;                               % Scaling in coeffs.
c_leg = scale.*(c_leg + c_rec);                        % Legendre coefficients.
c_leg = flipud(c_leg);
end

function C = constantOutTheFront(N) % (See Hale and Townsend, 2012)
%CONSTANTOUTTHEFRONT(N) = returns sqrt(4/pi)*gamma((0:N)+1)/gamma((0:N)+3/2))
% Initialise:
NN = (0:N)';
NN(1) = 1; % Set the first value different from 0 to avoid complications.
ds = -1/8./NN; s = ds; j = 1; ds(1) = 1;
while ( norm(ds(10:end)./s(10:end),inf) > eps/100 )
    j = j + 1;
    ds = -.5*(j-1)/(j+1)./NN.*ds;
    s = s + ds;
end
NN(1) = 0; % Reset the first value.
p2 = exp(s).*sqrt(4./(NN+.5)/pi);
% Stirling's series:
g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
    5246819/75246796800 -534703531/902961561600 ...
    -4483131259/86684309913600 432261921612371/514904800886784000];
eN = ones(N+1, 1); e9 = ones(1, 9);
ff1 = sum(bsxfun(@times, g, [eN, cumprod(bsxfun(@rdivide, e9, NN),2)]), 2);
ff2 = sum(bsxfun(@times, g, [eN, cumprod(bsxfun(@rdivide, e9, NN+.5),2)]), 2);
C = p2.*ff1./ff2;
% Use direct evaluation for the small values:
C(1:10) = sqrt(4/pi)*gamma((0:9)+1)./gamma((0:9)+3/2);
end

function v = dct1(c)
%DCT1   Compute a (scaled) DCT of type 1 using the FFT. 
% DCT1(C) returns T_N(X_N)*C, where X_N = cos(pi*(0:N))/N and T_N(X) = [T_0,
% T_1, ..., T_N](X) where T_k is the kth 1st-kind Chebyshev polynomial.
N = size(c, 1);                     % Number of terms.
ii = N-1:-1:2;                      % Indicies of interior coefficients.
c(ii) = 0.5*c(ii);                  % Scale interior coefficients.
v = ifft([c ; c(ii)]);              % Mirror coefficients and call FFT.
v = (N-1)*[ 2*v(N) ; v(ii) + v(2*N-ii) ; 2*v(1) ]; % Re-order.
v = flipud(v);                      % Flip the order.
end

function c = dst1Transpose(v, flag) %#ok<INUSD>
%DST1TRANSPOSE   Compute a transposed and scaled DST of type 1. 
% DST1TRANSPOSE(C) returns U_N(X)'*diag(sin(T_N))*C where T_N(k,1) = pi*(k-1)/N,
% k = 1:N+1, X_N = cos(T_N) and U_N(X) = [U_0, U_1, ..., U_N](X) where U_k is
% the kth 2nd-kind Chebyshev polynomial.
persistent Smat sint                % The same for each partition.
if ( nargin == 2 ), Smat = []; return, end % Clear persistent variables.
if ( isempty(Smat) ) % Construct conversion matrix:
    N = length(v) - 1;              % Degree of polynomial.
    dg = .5*ones(N-1, 1);           % Conversion matrix:
    Smat = spdiags([1 ; .5 ; dg], 0, N+1, N+1) + spdiags([0 ; 0 ; -dg], 2, N+1, N+1);
    sint = sin(pi*(0:N)'/N);        % Sin(theta).
end
c = (dct1(sint.*v)'/Smat)';         % Scaled DCT.
end
    
function w = cheb2_quadwts(n)       % Nick Hale, 2013 (See also Waldvogel, 2005)
c = 2./[1 , 1-(2:2:(n-1)).^2];      % Exact integrals of T_k (even)
c = [c, c(floor(n/2):-1:2)];        % Mirror for DCT via FFT 
w = ifft(c);                        % Interior weights
w([1,n]) = w(1)/2;                  % Boundary weights
end

function c_leg = cheb2leg_direct(c_cheb)
%CHEB2LEG_DIRECT   Convert Cheb to Leg coeffs using the 3-term recurrence.
N = length(c_cheb) - 1;             % Degree of polynomial.
if ( N <= 0 ), c_leg = c_cheb; return, end % Trivial case.
x = cos(.5*pi*(0:2*N)'/N);          % 2*N+1 Chebyshev grid (reversed order).
f = dct1([c_cheb ; zeros(N,1)]);    % Values on 2*N+1 Chebyshev grid.
w = cheb2_quadwts(2*N+1).';         % Clenshaw-Curtis quadrature weights.
Pm2 = 1; Pm1 = x;                   % Initialise.
L = zeros(2*N+1, N+1);              % Vandermonde matrix.
L(:,1:2) = [1+0*x, x];              % P_0 and P_1.
for k = 1:N-1 % Recurrence relation:
    P = ((2*k+1)/(k+1)*Pm1.*x-k/(k+1)*Pm2);           
    Pm2 = Pm1; Pm1 = P; 
    L(:,2+k) = P;
end
scale = (2*(0:N).'+1)/2;            % Scaling in coefficients.
c_leg = scale.*(L'*(w.*f));         % Legendre coefficients.
c_leg = flipud(c_leg);
end
