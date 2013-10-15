function c_cheb = leg2cheb(c_leg, M)
%LEG2CHEB   Convert Legendre coefficients to Chebyshev coefficients. 
%   C_CHEB = LEG2CHEB(C_LEG) converts the vector C_LEG of Legendre coefficients
%   to a vector C_CHEB of Chebyshev coefficients such that
%   C_CHEB(N)*T0 + ... + C_CHEB(1)*T{N-1} = C_LEG(N)*P0 + ... + C_LEG(1)*P{N-1}.

% Copyright 2013 Nick Hale and Alex Townsend, University of Oxford.

% This algorithm requires O( N(log N)^2 / log log N) operations and is based on
% rewritting an asymptotic formula for Legendre polynomials in a way that can be
% evaluated using discrete cosine transforms. For more details see:
%   N. Hale and A. Townsend, A fast, simple, and stable Chebyshev-Legendre
%   transform using an asymptotic formula, SISC (submitted) 2013.

c_leg = c_leg(:);                               % Make column vector.
c_leg = flipud(c_leg);                          % Lowest order coeffs first.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( nargin == 1 ), M = 10; end                 % No. of terms in expansion.
N = length(c_leg) - 1; NN = (0:N)';             % Degree of polynomial.
nM0 = min(floor(.5*(.25*eps*pi^1.5*gamma(M+1)/gamma(M+.5)^2)^(-1/(M+.5))),N);
aM = min(1/log(N/nM0), .5);                     % Block reduction factor (alpha)
K = ceil(log(N/nM0)/log(1/aM));                 % Number of block partitionss.

% Use direct approach if N is small:
if ( M == 0 || N < 513 || K == 0 ), c_cheb = leg2cheb_direct(c_leg); return, end

t = pi*(0:N)'/N;                                % Theta variable.   
nM = ceil(aM.^(K-1:-1:0)*N);                    % n_M for each block.
jK = zeros(K, 2);                               % Block locations in theta.
for k = 1:K
    tmp = find(t >= asin(nM0./nM(k)), 1) - 4;   % Where curve intersects aM^k*N.    
    jK(k,:) = [tmp+1, N+1-tmp];                 % Collect indicies.
end

%% %%%%%%%%%%%%%%%%%%% Recurrence / boundary region %%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_rec = zeros(N+1, 1);
jK2 = [jK(1,2)+1, N+1 ; jK(1:K-1,:) ; 1, N+1];
for k = 1:K % Loop over the partitions:
    j_bdy = [jK2(k+1,1):jK2(k,1)-1, jK2(k,2)+1:jK2(k+1,2)];
    x_bdy = cos(t(j_bdy));                      % Boundary indicies.
    tmp = c_leg(1) + c_leg(2)*x_bdy;            % Entries of mat-vec result.       
    Pm2 = 1; Pm1 = x_bdy;                       % Initialise recurrence.
    for n = 1:nM(k)-1  % Recurrence:
        P = (2-1/(n+1))*Pm1.*x_bdy-n/(n+1)*Pm2;
        Pm2 = Pm1; Pm1 = P;
        tmp = tmp + c_leg(n+2)*P;               % Update local LHS.
    end
    v_rec(j_bdy) = v_rec(j_bdy) + tmp;          % Global correction LHS.
end

%% %%%%%%%%%%%%%%%%% Asymptotics / interior region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_leg = constantOutTheFront(N).*c_leg;            % Scaling factor, eqn (3.3).
v_cheb = zeros(N+1, 1);                           % Initialise output vector.
dst1([], 1);                                      % Clear persistent storage.
for k = 1:K-1 % Loop over the block partitions:
    v_k = zeros(N+1, 1);                          % Initialise local LHS.
    hm = ones(N+1,1); hm([1:nM(k)+1,nM(k+1)+2:end]) = 0; % Initialise h_m.
    j_k = jK(k,1):jK(k,2);                        % t indicies of kth block.
    t_k = pi/2*ones(N+1,1); t_k(j_k) = t(j_k);    % Theta in kth block.
    denom = 1./sqrt(2*sin(t_k));                  % initialise denomenator.
    for m = 0:M-1 % Terms in asymptotic formula:
        denom = (2*sin(t_k)).*denom;              % Update denominator.
        u = sin((m+.5)*(.5*pi-t_k))./denom;       % Trig terms:
        v = cos((m+.5)*(.5*pi-t_k))./denom;
        hmc = hm.*c_leg;                          % h_M*c_leg.
        v_k = v_k + u.*dst1(hmc) + v.*dct1(hmc);  % Update using DCT1 and DST1.
        hm = ((m+0.5)^2./((m+1)*(NN+m+1.5))).*hm; % Update h_m.
    end
    v_cheb(j_k) = v_cheb(j_k) + v_k(j_k);         % Add terms to output vector.
end
dst1([], 1);                                      % Clear persistent storage.

%% %%%%%%%%%%%%%%%%%%%%%%%%% Combine for result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_cheb = v_cheb + v_rec;                          % Values on Chebyshev grid.
c_cheb = idct1(v_cheb);                           % Chebyshev coeffs.
c_cheb = flipud(c_cheb);

end

function C = constantOutTheFront(N) % (See Hale and Townsend, 2012)
%CONSTANTOUTTHEFRONT(N) returns sqrt(4/pi)*gamma((0:N)+1)/gamma((0:N)+3/2))
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

function v = dst1(c, flag) %#ok<INUSD>
%DST1   Discrete sine transform of type 1.
% DST1(C) returns diag(sin(T_N))*U_N(X)*C where T_N(k,1) = pi*(k-1)/N, k =
% 1:N+1, X_N = cos(T_N) and U_N(X) = [U_0, U_1, ..., U_N](X) where U_k is the
% kth 2nd-kind Chebyshev polynomial.
persistent Smat sint                % The same for each partition.
if ( nargin == 2 ), Smat = []; return, end % Clear persistent variables.
if ( isempty(Smat) ) % Construct conversion matrix:
    N = length(c) - 1;              % Degree of polynomial.
    dg = .5*ones(N-2, 1);           % Conversion matrix:
    Smat = spdiags([1 ; .5 ; dg], 0, N, N) + spdiags([0 ; 0 ; -dg], 2, N, N);
    sint = sin(pi*(0:N).'/N);       % Sin(theta).
end
v = sint.*dct1([Smat\c(2:end) ; 0]);% Scaled DCT.
end

function c = idct1(v)
%IDCT1   Convert values on a Cheb grid to Cheb coefficients (inverse DCT1).
% IDCT1(V) returns T_N(X_N)\V, where X_N = cos(pi*(0:N))/N and T_N(X) = [T_0,
% T_1, ..., T_N](X) where T_k is the kth 1st-kind Chebyshev polynomial.
N = size(v, 1);                        % Number of terms.
c = fft([v ; v(N-1:-1:2)])/(2*N-2);    % Laurent fold in columns and call FFT.
c = c(1:N);                            % Extract the first N terms.
if (N > 2), c(2:N-1) = 2*c(2:N-1); end % Scale interior coefficients.
if isreal(v), c = real(c); end         % Ensure a real output.
end

function c_cheb = leg2cheb_direct(c_leg)
%LEG2CHEB_DIRECT   Convert Leg to Cheb coeffs using the 3-term recurrence.
N = length(c_leg) - 1;              % Degree of polynomial.
if ( N <= 0 ), c_cheb = c_leg; return, end % Trivial case.
x = cos(pi*(0:N)'/N);               % Chebyshev grid (reversed order).
% Make the Legendre-Chebyshev Vandemonde matrix:
Pm2 = 1; Pm1 = x;                   % Initialise.
L = zeros(N+1, N+1);                % Vandermonde matrix. 
L(:,1:2) = [1+0*x, x];              % P_0 and P_1.     
for n = 1:N-1 % Recurrence relation:
    P = (2-1/(n+1))*Pm1.*x-n/(n+1)*Pm2;
    Pm2 = Pm1; Pm1 = P; 
    L(:,2+n) = P;
end
v_cheb = L*c_leg;                   % Values on Chebyshev grid.
c_cheb = idct1(v_cheb);             % Chebyshev coefficients.
c_cheb = flipud(c_cheb);
end
