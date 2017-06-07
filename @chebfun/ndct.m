function vals = ndct(x, coeffs, th)
%NDCT   Nonuniform discrete Chebyshev transform
%   V = NDCT(X, C) computes the nonuniform Chebyshev transform, which evaluates
%   the Chebyshev expansion with coefficients C at the points X, i.e.,
%
%             V_j = sum_k  C(k) T_k( X(j) ).
%
%   X must be a column vector, but C may be a column vector or a matrix.
%
%   V = NDCT(~, C, TH) is similar, but with X(j) in the expression above given
%   by X(j) = cos(T(j)) so that 
%
%             V_j = sum_k  C(k) T_k( cos(TH(j) ).
%
%   This is useful as TH(j) can often be computed more accurately than X(j).
%
%   V = NDCT(C) is the same as evaluating the Chebyshev series at Legendre 
%   nodes, i.e., 
%     [~, ~, ~, TH] = legpts( size(C,1) );  
%     V = NDCT(~, C, TH);
%
% See also CHEBTECH.clenshaw, CHEBTECH2.coeffs2vals, CHEBTECH1.coeffs2vals.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: The algorithm in this MATLAB script is based on the paper
%   [1] D. Ruiz--Antolin and A. Townsend, "A nonuniform fast Fourier transform
%   based on low rank approximation", submitted, 2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For backwards compatibility with chebfun.ndct(c):
if ( nargin == 1 ) 
    coeffs = x; 
    [x,~,~,th] = legpts( size(coeffs, 1) );
    vals = chebfun.ndct(x, coeffs, th);
    return
end

% n is the length of the chebtech. If m > 1 then it is array-valued.
[n, m] = size(coeffs);

% Convert evaluation points to theta space:
if ( nargin < 3 )
    % Compute theta values of evaluation points:  
    th = real(acos(x)/2/pi);
else
    % Scale to [0 1];
    th = real(th/2/pi);
end

% th should be a column vector.
if ( size(th, 2) > 1 )
    warning('CHEBFUN:CHEBTECH:fastChebyshevEval:xDim', ...
        'Evaluation points should be a column vector.');
    th = th(:);
end

% Convert a NUDCT into a NUFFT by mirroring (see chebtech2.coeffs2vals):
c = coeffs;
c(2:end,:) = c(2:end,:)/2;
c = [c(end:-1:2,:) ; c];
% Get new length of vector: 
N = 2*n - 1;

% We always take K = 16 terms in the low-rank approximation.
K = 16; 

% Closest points: 
s = round(N*th);
t = mod(s, N) + 1;
ds = 2*(N*th - s);

% Construct a low rank approximation NUDFT./DFT: 
U = ChebP(K-1, ds) * besselCoeffs(K);
k = 2*(-(n-1):(n-1))'/(2*n-1);
V = ChebP(K-1, k);

vals = zeros(numel(x), m);
% Loop over the columns of c for array-valued input:
for l = 1:m
    % The NUDFT can now be written as a sum of diagonally-scaled DFTs:
    C = repmat(c(:,l), 1, K);
    tmp1 = ifftshift(C.*V, 1);
    tmp2 = fft(tmp1, [], 1);
    vals(:,l) = sum(U.*tmp2(t,:), 2);
end

% If the coefficients were real/imaginary, then enforce it on vals: 
%  (Note: We don't bother to check each column of coeffs for real/imag.)
if ( isreal(coeffs) )
    vals = real(vals);
elseif ( isreal(1i*coeffs) )
    vals = imag(vals);
end

end

function T = ChebP(n, x)
% Evaluate Chebyshev polynomials of degree 0,...,n at points in x. Use the
% three-term recurrence relation:
N = size(x, 1);
T = zeros(N, n+1);
T(:,1) = 1;
if ( n == 0 )
    return
end
T(:,2) = x;
twoX = 2*x;
for k = 2:n
    T(:,k+1) = twoX.*T(:,k) - T(:,k-1);
end
end

function cfs = besselCoeffs(K)
% Bivarate Chebyshev coefficients for the function f(x,y) = exp(-i*x.*y) on the
% domain [-1/2, 1/2]x[0,2*pi] are given by Lemma A.2 of Townsend's DPhil thesis.

% Here, we tabulate them so that they are extremely fast to compute!: 
cfs = [ 0.725276916440514 0 -0.263810811846140 0 0.010721845410224 0 -0.000188568764214 0 0.000001845983729 0 -0.000000011505371 0 0.000000000049650 0 -0.000000000000157 0 ;
      0 -1.237209415222620i 0 0.106368016662154i 0 -0.002843803888523i 0 0.000037314601460i 0 -0.000000291470262i 0 0.000000001511615i 0 -0.000000000005586i 0 0.000000000000015i;
      -0.263810811846140 0 -0.249420239398158 0 0.014106236744931 0 -0.000281370589589 0 0.000002945880969 0 -0.000000019147182 0 0.000000000085036 0 -0.000000000000275 0;
      0 0.106368016662154i 0 0.033077433013562i 0 -0.001395694044102i 0 0.000022213402600i 0 -0.000000193519978i 0 0.000000001077127i 0 -0.000000000004183i 0 0.000000000000012i;
      0.010721845410224 0 0.014106236744931 0 0.003272735109015 0 -0.000110186049486 0 0.000001459236547 0 -0.000000010886489 0 0.000000000052987 0 -0.000000000000183 0 ; 
      0 -0.002843803888523i 0 -0.001395694044102i 0 -0.000258373068366i 0 0.000007238310730i 0 -0.000000082089523i 0 0.000000000535538i 0 -0.000000000002316i 0 0.000000000000007i;
      -1.885687642135953e-04 0 -2.813705895892196e-04 0 -1.101860494855995e-04 0 -1.697297037000288e-05 0 4.071920170840575e-07 0 -4.038224124439603e-09 0 2.340743035825754e-11 0 -9.107573296934322e-14 0; 
      0 3.731460145969448e-05i 0 2.221340260014939e-05i 0 7.238310730031504e-06i 0 9.548164341984993e-07i 0 -2.003096829508016e-08i 0 1.765036265437698e-10i 0 -9.204996298527833e-13i 0 3.255202808838026e-15i;
      1.845983728936490e-06 0 2.945880968932915e-06 0 1.459236547111028e-06 0 4.071920170840575e-07 0 4.697021778082508e-08 0 -8.755181580605211e-10 0 6.941023444886516e-12 0 -3.290023459530936e-14 0 ;
      0 -2.914702622193488e-07i 0 -1.935199775867226e-07i 0 -8.208952270529682e-08i 0 -2.003096829508016e-08i 0 -2.052985055408922e-09i 0 3.442984249400081e-11i 0 -2.480840754980324e-13i 0 1.077723939077946e-15i;
      -1.150537142155097e-08 0 -1.914718184703716e-08 0 -1.088648898321836e-08 0 -4.038224124439603e-09 0 -8.755181580605211e-10 0 -8.073385051984342e-11 0 1.230581586777325e-12 0 -8.126572662991549e-15 0 ; 
      0 1.511615196341062e-09i 0 1.077126955246971e-09i 0 5.355382878799667e-10i 0 1.765036265437698e-10i 0 3.442984249400081e-11i 0 2.885566203117638e-12i 0 -4.031057077166038e-14i 0 2.456927663965262e-16i;
      4.965029850164547e-11 0 8.503608974664060e-11 0 5.298703065162090e-11 0 2.340743035825754e-11 0 6.941023444886516e-12 0 1.230581586777325e-12 0 9.452345289165531e-14 0 -1.218719878432285e-15 0; 
      0 -5.586166703739336e-12i 0 -4.183174389936338e-12i 0 -2.315969292837316e-12i 0 -9.204996298527833e-13i 0 -2.480840754980324e-13i 0 -4.031057077166038e-14i 0 2.857751919953105e-15i 0 0 ; 
      -1.571252307825089e-13 0 -2.747999062823869e-13 0 -1.828391460048661e-13 0 -9.107573296934322e-14 0 -3.290023459530936e-14 0 -8.126572662991549e-15 0 -1.218719878432285e-15 0 0 0;
      0 1.545890088268538e-14i 0 1.201101735269839e-14i 0 7.190168405679147e-15i 0  3.255202808838026e-15i 0 1.077723939077946e-15i 0 2.456927663965262e-16i 0 0 0 0];  
cfs = cfs(1:K,1:K);
end
