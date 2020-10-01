function vals = fevalm( f, r, lam, th)
%FEVALM   Evaluate a BALLFUN in spherical coordinates.
% 
% Z = FEVALM(F, R, LAM, TH) returns a matrix of values Z of size
% length(R)-by-length(TH)-by-length(LAM). (R,LAM,TH) are the spherical 
% coordinates for the evaluation points in the ball, with 0<=R<=1, 
% -pi <= LAM <= pi the azimuthal angle and 0 <= TH <= pi the elevation 
% (polar) angle (both measured in radians). R, LAM, and TH should be 
% vectors of doubles. This is equivalent to making a ndgrid of the vectors 
% R, LAM, and TH and then using FEVAL to evaluate at that grid.
%
% See also FEVAL. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    vals = [];
    return
end

F = f.coeffs;
[m, n, p] = size( f );

% Get the size of the lists
Nr = length( r );
Nlam = length( lam );
Nth = length( th );

% Transform the lists to vectors
r = reshape(r, Nr, 1);
lam = reshape(lam, Nlam, 1);
th = reshape(th, Nth, 1);

% Evaluate f at the points r
F1 = reshape(F, m, n*p);
G = reshape(clenshaw_vec( r, F1 ), Nr, n, p);
G = permute(G, [2, 1, 3]);
G = reshape(G, n, Nr*p);

% Permute G to evaluate f at lambda
H = reshape(horner_vec_cmplx( lam/pi, G ), Nlam, Nr, p);
H = permute(H, [3, 2, 1]);
H = reshape(H, p, Nr*Nlam);

% Evaluate H at theta
vals = reshape(horner_vec_cmplx( th/pi, H ), Nth, Nr, Nlam); 

% Permute H to get the array of values r x lambda x theta
vals = permute(vals, [2, 3, 1]);

% Return real values if the function is real
if f.isReal
   vals = real(vals); 
end
end

% Pulled from chebtech/clenshaw: 
function y = clenshaw_vec(x, c)
% Clenshaw scheme for array-valued functions.
x = repmat(x(:), 1, size(c, 2));
bk1 = zeros(size(x, 1), size(c, 2)); 
bk2 = bk1;
e = ones(size(x, 1), 1);
x = 2*x;
n = size(c, 1)-1;
for k = (n+1):-2:3
    bk2 = e*c(k,:) + x.*bk1 - bk2;
    bk1 = e*c(k-1,:) + x.*bk2 - bk1;
end
if ( mod(n, 2) )
    [bk1, bk2] = deal(e*c(2,:) + x.*bk1 - bk2, bk1);
end
y = e*c(1,:) + .5*x.*bk1 - bk2;
end

% Pulled from trigtech/horner: 
% Vectorized version of Horner's scheme for evaluating multiple
% polynomials of the same degree at the same locations.
function q = horner_vec_cmplx(x, c)
nValsX = size(x, 1);
N = size(c, 1);
% Just return the constant term.
if N == 1
    q = repmat(c(N,:),[nValsX 1]);
    return
end

z = exp(1i*pi*x); % Use complex exponential
numCols = size(c,2);
if nValsX > 1
    z = repmat(z,[1 numCols]);
    e = ones(nValsX,1);
    q = repmat(c(N,:),[nValsX 1]);
    for j = N-1:-1:2
        q = e*c(j,:) + z.*q;
    end
    if mod(N, 2) == 1
        q = bsxfun(@times, exp(-1i*pi*(N-1)/2*x), e*c(1,:) + z.*q);
    else
        q = bsxfun(@times, exp(-1i*pi*(N/2-1)*x), q) + cos(N/2*pi*x)*c(1,:);
    end
else
    q = c(N,:);
    for j = N-1:-1:2
        q = z*q + c(j,:);
    end
    if mod(N, 2) == 1
        q = exp(-1i*pi*(N-1)/2*x)*(z.*q + c(1,:));
    else
        q = exp(-1i*pi*(N/2-1)*x)*q + cos(N/2*pi*x)*c(1,:);
    end
end

end