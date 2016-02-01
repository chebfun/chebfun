function c_cheb = leg2cheb(c_leg, varargin)
%LEG2CHEB   Convert Legendre coefficients to Chebyshev coefficients. 
%   C_CHEB = LEG2CHEB(C_LEG) converts the vector C_LEG of Legendre coefficients
%   to a vector C_CHEB of Chebyshev coefficients such that 
%       C_CHEB(1)*T0 + ... + C_CHEB(N)*T{N-1} = ...
%           C_LEG(N)*P0 + ... + C_LEG(1)*P{N-1}, 
%   where P{k} is the degree k Legendre polynomial normalized so that max(|P{k}|
%   = 1.
% 
%   C_CHEB = LEG2CHEB(C_LEG, 'norm') is as above, but with the Legendre
%   polynomials normalized to be orthonormal.
%
%   C = LEG2CHEB(C_LEG, 'trans') returns the `transpose' of the LEG2CHEB
%   operator applied to C_LEG. That is, if C_CHEB = B*C_LEG, then C = B'*C_LEG.
%
%   If C_LEG is a matrix then the LEG2CHEB operation is applied to each column.
%
% See also CHEB2LEG. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  This algorithm requires O( N(log N)^2 / log log N) operations and is based on
%  rewritting an asymptotic formula for Legendre polynomials in a way that can
%  be evaluated using discrete cosine transforms. For more details see:
%   N. Hale and A. Townsend, A fast, simple, and stable Chebyshev-Legendre
%   transform using an asymptotic formula, SISC, 36 (2014), pp. A148-A167.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 7;                                          % No. of terms in expansion.
normalize = 0;                                  % Default - no normalize.
trans = 0;                                      % Default - no transpose.
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
[N, n] = size(c_leg);                           % Number of columns.
% Do normalization:
if ( normalize ) 
    c_leg = bsxfun(@times, c_leg, sqrt((0:N-1)'+1/2) ); 
end
% Trivial case:
if ( N < 2 )
    c_cheb = c_leg;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = N - 1; NN = (0:N)';                         % Degree of polynomial.
nM0 = min(floor(.5*(.25*eps*pi^1.5*gamma(M+1)/gamma(M+.5)^2)^(-1/(M+.5))), N);
aM = min(1/log(N/nM0), .5);                     % Block reduction factor (alpha)
K = ceil(log(N/nM0)/log(1/aM));                 % Number of block partitions.
 
if ( M == 0 || N < 513 || K == 0 )

    %%%%%%%%%%%%%%%%%%%% Use direct approach if N is small %%%%%%%%%%%%%%%%%%%%%
    L = legchebvandermonde(N);
    if ( ~trans )
        c_cheb = idct1(L*c_leg);
    else
        c_cheb = L.'*idct1(c_leg);
    end
    
else

    %%%%%%%%%%%%%%%%%%%% Use asymptotics for large N %%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = pi*NN/N;                                    % Theta variable.   
    nM = ceil(aM.^(K-1:-1:0)*N);                    % n_M for each block.
    jK = zeros(K, 2);                               % Block locations in theta.
    for j = 1:K
        tmp = find(t >= asin(nM0./nM(j)), 1) - 4;   % Curve intersects aM^k*N.    
        jK(j,:) = [tmp+1, N+1-tmp];                 % Collect indicies.
    end
    jK2 = [jK(1,2)+1, N+1 ; jK(1:K-1,:) ; 1, N+1];

    if ( ~trans )
        c_cheb = leg2cheb_fast();
    else
        c_cheb = leg2cheb_transpose_fast();
    end
    
    if ( isreal(c_leg) )
        c_cheb = real(c_cheb);
    elseif ( isreal(1i*c_leg) )
        c_cheb = imag(c_cheb);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FAST METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regular:

function c_cheb = leg2cheb_fast()
%%%%%%%%%%%%%%%%%%%%%% Recurrence / boundary region %%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_rec = zeros(N+1, n);
for k = 1:K % Loop over the partitions:
    j_bdy = [jK2(k+1,1):jK2(k,1)-1, jK2(k,2)+1:jK2(k+1,2)];
    x_bdy = cos(t(j_bdy));                      % Boundary indicies.
    P = zeros(length(x_bdy), nM(k)+1); P(:,1) = 1; P(:,2) = x_bdy; % Initialise.
    for kk = 1:(nM(k)-1)                        % Recurrence:
        P(:,kk+2) = (2-1/(kk+1))*P(:,kk+1).*x_bdy - (1-1/(kk+1))*P(:,kk);
    end
    tmp = P*c_leg(1:nM(k)+1,:);
    v_rec(j_bdy,:) = v_rec(j_bdy,:) + tmp;      % Global correction LHS.
end
%%%%%%%%%%%%%%%%%%%% Asymptotics / interior region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = constantOutTheFront(N);
c_leg = bsxfun(@times, c_leg, C);                 % Scaling factor, eqn (3.3).
v_cheb = zeros(N+1, n);                           % Initialise output vector.
for k = 1:K-1 % Loop over the block partitions:
    v_k = zeros(N+1, n);                          % Initialise local LHS.
    hm = ones(N+1,n); hm([1:nM(k)+1, nM(k+1)+2:end],:) = 0; % Initialise h_m.
    j_k = jK(k,1):jK(k,2);                        % t indicies of kth block.
    t_k = pi/2*ones(N+1, 1); t_k(j_k) = t(j_k);   % Theta in kth block.
    denom = 1./sqrt(2*sin(t_k));                  % initialise denomenator.
    for m = 0:(M-1)                               % Terms in asymptotic formula.
        denom = (2*sin(t_k)).*denom;              % Update denominator.
        u = sin((m+.5)*(.5*pi-t_k))./denom;       % Trig terms:
        v = cos((m+.5)*(.5*pi-t_k))./denom;
        hmc = bsxfun(@times, c_leg,hm);           % h_M*c_leg.
        % Update using DCT1 and DST1:
        U = dst1(hmc(2:end,:));
        Uu = bsxfun(@times, U, u);
        V = dct1(hmc);
        Vv = bsxfun(@times, V, v);
        v_k = v_k + Uu + Vv; 
        hm = bsxfun(@times, hm, ((m+0.5)^2./((m+1)*(NN+m+1.5)))); % Update h_m.
    end
    v_cheb(j_k,:) = v_cheb(j_k,:) + v_k(j_k,:);   % Add terms to output vector.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Combine for result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_cheb = v_cheb + v_rec;                          % Values on Chebyshev grid.
c_cheb = idct1(v_cheb);                           % Chebyshev coeffs.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transpose:

function c = leg2cheb_transpose_fast()
f = idct1(c_leg);
%%%%%%%%%%%%%%%%%%%%%% Recurrence / boundary region %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_rec = zeros(N+1, n);
for k = 1:K    % Loop over the block partitions:
    j_bdy = [jK2(k+1,1):jK2(k,1)-1, jK2(k,2)+1:jK2(k+1,2)];% Boundary indicies.
    x_bdy = cos(t(j_bdy));                         % Boundary x values.
    f_bdy = f(j_bdy,:);                            % f at boundary.
    P = zeros(size(f_bdy,1), nM(k)-1); 
    P(:,1) = 1; P(:,2) = x_bdy; % Initialise.
    for kk = 1:(nM(k)-1)        % Recurrence:
        P(:,kk+2) = (2-1/(kk+1))*(P(:,kk+1).*x_bdy) - (1-1/(kk+1))*P(:,kk); 
    end
    ck = P'*f_bdy;
    idxk = 1:(nM(k)+1);
    c_rec(idxk,:) = c_rec(idxk,:) + ck;            % Global correction LHS.
end
%%%%%%%%%%%%%%%%%%%% Asymptotics / interior region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = constantOutTheFront(N);                   % Scaling in asymptotic expansion.
c = zeros(N+1, n);                                 % Initialise output.
for k = 1:K-1 % Loop over the block partitions.
    c_k = zeros(N+1, n);                           % Initialise local LHS.
    hm = ones(N+1,n); hm([1:nM(k)+1,nM(k+1)+2:end],:) = 0; % Initialise h_m.
    j_k = jK(k,1):jK(k,2);                         % t indices of kth block.
    t_k = pi/2 + 0*t; t_k(j_k) = t(j_k);           % t indices of kth block.
    f_k = 0*f; f_k(j_k,:) = f(j_k,:);              % RHS f of kth block.
    denom = 1./sqrt(2*sin(t_k));                   % Initialise denominator.
    for m = 0:(M-1)                                % Terms in asymptotic formula
        denom = (2*sin(t_k)).*denom;               % Update denominator.
        u = sin((m+.5)*(.5*pi-t_k))./denom;        % Trig terms:
        v = cos((m+.5)*(.5*pi-t_k))./denom;
        Dv = bsxfun(@times, f_k, v);
        Tv = dct1(Dv);                             % Compute T'v.
        Du = bsxfun(@times, f_k, u);
        Uu = dst1Transpose(Du);                    % Compute U'u*sin(t).
        c_k = c_k + hm.*([zeros(1,n);Uu(1:N,:)] + Tv(1:N+1,:));   % Update LHS.
        hm = bsxfun(@times, hm, ((m+0.5)^2./((m+1)*(NN+m+1.5)))); % Update h_m.
    end
    c = c + bsxfun(@times, c_k, C);                % Append to global LHS.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Combine for result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = c + c_rec;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% DIRECT METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = legchebvandermonde(N)
x = cos(pi*(0:N)'/N);                       % Chebyshev grid (reversed order).
% Make the Legendre-Chebyshev Vandemonde matrix:
Pm2 = 1; Pm1 = x;                           % Initialise.
L = zeros(N+1, N+1);                        % Vandermonde matrix. 
L(:,1:2) = [1+0*x, x];                      % P_0 and P_1.     
for n = 1:N-1                               % Recurrence relation:
    P = (2-1/(n+1))*Pm1.*x - (1-1/(n+1))*Pm2;  
    Pm2 = Pm1; Pm1 = P; 
    L(:,2+n) = P;
end
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

function v = dst1(c)
%DST1   Discrete sine transform of type 1.
% DST1(C) returns diag(sin(T))*U(cos(T))*C where T(k) = pi*k/(N+1), k = 0:N+1,
% and U(x) = [U_0, U_1, ..., U_N}](x) (where U_k is the kth 2nd-kind Chebyshev
% polynomial), and N = length(C) - 1.

z = zeros(1, size(c, 2));                   % Padding;
v = [z ; chebfun.dst(c(1:end-1,:), 1) ; z]; % DST-I.

end

function c = idct1(v)
%IDCT1   Convert values on a Cheb grid to Cheb coefficients (inverse DCT1).
% IDCT1(V) returns T(X)\V, where X = cos(pi*(0:N)/N), T(X) = [T_0, T_1, ...,
% T_N](X) (where T_k is the kth 1st-kind Chebyshev polynomial), and N =
% length(V) - 1.

c = chebfun.idct(v, 1);                     % IDCT-I
c([1,end],:) = .5*c([1,end],:);             % Scale.

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
%CONSTANTOUTTHEFRONT(N) returns sqrt(4/pi)*gamma((0:N)+1)/gamma((0:N)+3/2))
NN = (0:N).';
C = sqrt(4/pi)*exp(gammaln(NN+1) - gammaln(NN+3/2));
end
