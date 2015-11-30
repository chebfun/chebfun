function c_leg = cheb2leg(c_cheb, varargin)
%CHEB2LEG   Convert Chebyshev coefficients to Legendre coefficients. 
%   C_LEG = CHEB2LEG(C_CHEB) converts the vector C_CHEB of Chebyshev
%   coefficients to a vector C_LEG of Legendre coefficients such that
%       C_CHEB(1)*T0 + ... + C_CHEB(N)*T{N-1} = ...
%           C_LEG(1)*P0 + ... + C_LEG(N)*P{N-1},
%   where P{k} is the degree k Legendre polynomial normalized so that
%   max(|P{k}|) = 1.
% 
%   C_LEG = CHEB2LEG(C_CHEB, 'norm') is as above, but with the Legendre
%   polynomials normalized to be orthonormal.
%
%   If C_CHEB is a matrix then the CHEB2LEG operation is applied to each column.
%
% See also LEG2CHEB.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
% This algorithm requires O( N(log N)^2 / log log N) operations and is based on
% rewritting an asymptotic formula for Legendre polynomials in a way that can be
% evaluated using discrete cosine transforms. For more details see:
%   N. Hale and A. Townsend, A fast, simple, and stable Chebyshev-Legendre
%   transform using an asymptotic formula, SISC, 36 (2014), pp. A148-A167.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 7;                                         % No. of terms in expansion.
normalize = 0;                                 % Default - no normalize.
trans = 0;                                     % Default - no transpose.
for j = 1:numel(varargin)
    if ( strncmpi(varargin{j}, 'norm', 4) )
        normalize = 1;
    elseif ( strncmpi(varargin{j}, 'trans', 4) )
        trans = 1;
    end
end
if ( normalize && trans )
    error('CHEBFUN:CHEBFUN:LEG2CHEB:normtrans', ...
        'No support for both ''norm'' and ''trans'' in LEG2CHEB.')
end

[N, n] = size(c_cheb);                         % Number of columns.
% Trivial case:
if ( N < 2 )
    c_leg = c_cheb;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = N - 1; NN = (0:N).';                      % Degree of polynomial.
nM0 = min(floor(.5*(.25*eps*pi^1.5*gamma(M+1)/gamma(M+.5)^2)^(-1/(M+.5))), N);
aM = min(1/log(N/nM0), .5);                   % Block reduction factor (alpha_M)
K = ceil(log(N/nM0)/log(1/aM));               % Number of block partitions

% Use direct approach if N is small:
if ( M == 0 || N < 513 || K == 0 ) 
    c_leg = cheb2leg_direct(c_cheb); 
    if ( normalize ), 
        c_leg  = bsxfun(@times, c_leg, 1./sqrt(NN+.5) ); 
    end
    return 
end

f = dct1([c_cheb ; zeros(N,n)]);              % Values on a 2*N+1 Cheb grid.
w = chebtech2.quadwts(2*N+1);                 % C-C quadrature weights.
wf = bsxfun(@times, f, w.');                  % Scale f by C-C weights.
t = pi*(0:2*N)'/(2*N);                        % 2*N+1 theta grid.
nM = ceil(aM.^(K-1:-1:0)*N);                  % n_M for each block.
jK = zeros(K, 2);                             % Block locations in theta.
for k = 1:K    % Find where curve intersects a^k*N:
    tmp = find(t >= asin(nM0./nM(k)), 1) - 4; % Where curve intersects aM^k*N.    
    jK(k,:) = [tmp+1, 2*N+1-tmp];             % Collect indicies.
end
C = constantOutTheFront(N);                   % Scaling in asymptotic expansion.
nM(end) = N+2; % For convenience (avoids treating final block differently).

%%%%%%%%%%%%%%%%%%%%%% Recurrence / boundary region %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_rec = zeros(N+1, n);
jK2 = [jK(1,2)+1, 2*N+1 ; jK(1:K-1,:) ; 1, 2*N+1];
for k = 1:K    % Loop over the block partitions:
    j_bdy = [jK2(k+1,1):jK2(k,1)-1, jK2(k,2)+1:jK2(k+1,2)];% Boundary indicies.
    x_bdy = cos(t(j_bdy));                             % Boundary x values.
    wf_bdy = wf(j_bdy,:);                              % w.*f at boundary.
    tmp = zeros(N+1,n);                                % Local LHS storage.
    tmp(1:2,:) = [sum(wf_bdy) ; x_bdy'*wf_bdy];        % First two terms.
    Pm2 = 1; Pm1 = x_bdy;                              % Initialise recurrence.
    for kk = 1:min(nM(k)-3, N-1) % Recurrence:
        P = (2-1/(kk+1))*Pm1.*x_bdy - (1-1/(kk+1))*Pm2; 
        Pm2 = Pm1; Pm1 = P;
        tmp(kk+2,:) = P'*wf_bdy;                       % Update local LHS.
    end
    c_rec = c_rec + tmp;                               % Global correction LHS.
end
%%%%%%%%%%%%%%%%%%%% Asymptotics / interior region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_leg = zeros(N+1, n);                                 % Initialise output.
for k = 1:K-1 % Loop over the block partitions.
    c_k = zeros(N+1, n);                               % Initialise local LHS.
    hm = ones(N+1,n); hm([1:nM(k)-1,nM(k+1):end],:) = 0; % Initialise h_m.
    j_k = jK(k,1):jK(k,2);                             % t indices of kth block.
    t_k = pi/2 + 0*t; t_k(j_k) = t(j_k);               % t indices of kth block.
    wf_k = 0*wf; wf_k(j_k,:) = wf(j_k,:);              % RHS w.*f of kth block.
    denom = 1./sqrt(2*sin(t_k));                       % Initialise denominator.
    for m = 0:M-1 % Terms in the asymptotic formula:
        denom = (2*sin(t_k)).*denom;                   % Update denominator.
        u = sin((m+.5)*(.5*pi-t_k))./denom;            % Trig terms:
        v = cos((m+.5)*(.5*pi-t_k))./denom;
        Dv = bsxfun(@times, wf_k, v);
        Tv = dct1(Dv);                                 % Compute T'v.
        Du = bsxfun(@times, wf_k, u);
        Uu = dst1Transpose(Du);                        % Compute U'u*sin(t).
        c_k = c_k + hm.*([zeros(1,n);Uu(1:N,:)] + Tv(1:N+1,:));% Update LHS.
        hm = bsxfun(@times, hm, ((m+0.5)^2./((m+1)*(NN+m+1.5))));  % Update h_m.
    end
    c_leg = c_leg + bsxfun(@times, c_k, C);            % Append to global LHS.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Combine for result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_leg = c_leg + c_rec;
scl = (NN+.5);                                       % Scaling in coeffs.
if ( normalize )
    c_leg  = bsxfun(@times, c_leg, sqrt(scl) ); 
else
    c_leg  = bsxfun(@times, c_leg, scl); 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% DIRECT METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c_leg = cheb2leg_direct(c_cheb)
%CHEB2LEG_DIRECT   Convert Cheb to Leg coeffs using the 3-term recurrence.
[N, m] = size(c_cheb);              % Number of columns.
N = N - 1;                          % Degree of polynomial.
if ( N <= 0 ), c_leg = c_cheb; return, end % Trivial case.
x = cos(.5*pi*(0:2*N).'/N);         % 2*N+1 Chebyshev grid (reversed order).
f = dct1([c_cheb ; zeros(N, m)]);   % Values on 2*N+1 Chebyshev grid.
w = chebtech2.quadwts(2*N+1).';     % Clenshaw-Curtis quadrature weights.
Pm2 = 1; Pm1 = x;                   % Initialise.
L = zeros(2*N+1, N+1);              % Vandermonde matrix.
L(:,1:2) = [1+0*x, x];              % P_0 and P_1.
for k = 1:N-1 % Recurrence relation:
    P = (2-1/(k+1))*Pm1.*x - (1-1/(k+1))*Pm2;    
    Pm2 = Pm1; Pm1 = P; 
    L(:,2+k) = P;
end
scale = (2*(0:N).'+1)/2;            % Scaling in coefficients [NIST, (18.3.1)].
c_leg = bsxfun(@times, L.'*(bsxfun(@times, f ,w)), scale); % Legendre coeffs.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% DCT METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = dct1(c)
%DCT1   Compute a (scaled) DCT of type 1 using the FFT. 
% DCT1(C) returns T(X)*C, where X = cos(pi*(0:N)/N) and T(X) = [T_0, T_1, ...,
% T_N](X) (where T_k is the kth 1st-kind Chebyshev polynomial), and N =
% length(C) - 1;

c([1,end],:) = 2*c([1,end],:);              % Scale.
v = chebfun.dct(c, 1);                      % DCT-I.

end

function c = dst1Transpose(v)
%DST1TRANSPOSE   Compute a transposed and scaled DST of type 1. 
% DST1TRANSPOSE(C) returns U(cos(T))'*diag(sin(T))*V where T(k,1) = pi*(k-1)/N,
% k = 1:N+1, Xand U_N(X) = [U_0, U_1, ..., U_N](X) (where U_k is the kth
% 2nd-kind Chebyshev polynomial), and N = legnth(V) - 1.

m = size( v, 2 ); 
c = [ chebfun.dst( v(2:end-1, :), 1 ) ; zeros(2, m) ]; 
c(end,:) = -c(end-2,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = constantOutTheFront(N) % (See Hale and Townsend, 2014)
%CONSTANTOUTTHEFRONT(N) = returns sqrt(4/pi)*gamma((0:N)+1)/gamma((0:N)+3/2))
NN = (0:N).';
C = sqrt(4/pi)*exp(gammaln(NN+1) - gammaln(NN+3/2));
end
