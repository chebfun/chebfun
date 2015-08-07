function f = diff( f, varargin )
% DIFF    Derivative of a spherefun in Cartesian coordinates.
%
%  F = DIFF( F ) computes the first derivative of F with respect to x.
%
%  F = DIFF( F, DIM )  computes the first derivative of F. If DIM = 1, the
%  derivative is taken in the x-direction. If DIM = 2, the derivative
%  is taken in the y-direction and if DIM = 3, the derivative is taken in
%  the z-direction.
%
%  F = DIFF( F, DIM, K) computes the kth derivatives of F in the variable
%  given by DIM.
%
%  See also GRADIENT, LAPLACIAN

% Parse user inputs:
if ( nargin == 1 )
    dim = 1;
    K = 1;
elseif ( nargin == 2 )
    dim = varargin{1};
    K = 1;
else
    dim = varargin{1};
    K = varargin{2};
end

if K > 1
    error('SPHEREFUN:DIFF:ORDER', 'Only first derivatives currently allowed.');
end

if ( dim ~= 1 && dim ~= 2 && dim ~= 3 )
    error('SPHEREFUN:DIFF:DIM', 'Unrecognized coordinate dimension');
end

if ( abs( K - round(K) ) > eps )
    error('SPHEREFUN:DIFF:DIFFORDER', 'Fractional derivatives not allowed')
end
K = round( K );

% Track whether f is real.
realf = isreal(f);

% We are going to work at the tech level to make things faster.
[cols, D, rows] = cdr( f );

coltechs = cols.funs{1}.onefun;
rowtechs = rows.funs{1}.onefun;

% Compute the derivatives
d_coltechs = diff(coltechs)/pi;
d_rowtechs = diff(rowtechs)/pi;

% We will do everyting in value space on [-1,1] at the tech level.

% Evaluate at the half grid points in theta, so the poles are not included.
m = length(cols);
n = length(rows);

h = 2*pi/m;
shift = h/2/pi;

coltechs = circshift(coltechs,-shift);
colv = coltechs.values;
d_coltechs = circshift(d_coltechs,-shift);
dcolv = d_coltechs.values;

% Theta at half-grid points (on on [-1,1]);
th = trigpts(m,[-1,1]) + h/2/pi;
lam = trigpts(n,[-1,1]).';

% Evalute rows and the derivatives at the grid
rowv = rowtechs.values.';
drowv = d_rowtechs.values.';

if ( dim == 1 )            % x
    val = -((1./sin(pi*th))*sin(pi*lam)).*(colv*D*drowv) + ...
        (cos(pi*th)*cos(pi*lam)).*(dcolv*D*rowv);
elseif ( dim == 2)         % y
    val = ((1./sin(pi*th))*cos(pi*lam)).*(colv*D*drowv) + ...
        (cos(pi*th)*sin(pi*lam)).*(dcolv*D*rowv);
else
    val = -bsxfun(@times,dcolv*D*rowv,sin(pi*th));
end

%
% Shift back to regular grid points in theta.
%

% Just do the works ourselves rather than call trigtech as the code will be
% faster and we know how to easily shift back by h/2.

n = length(rows);
idcol1 = m:-1:m/2+1;
idcol2 = m/2:-1:1;
% Enforce exact symmetry
val = 0.5*(val + [val(idcol1,n/2+1:n) val(idcol1,1:n/2) ; ...
                  val(idcol2,n/2+1:n) val(idcol2,1:n/2)] );

% Wave numbers ordering according to MATLAB's FFT
m1 =  floor((m-1)/2);
m2 = (m/2)*ones(rem(m+1,2));
waveNum = [(0:m1)  m2 (-m1:-1)]';

val = ifft(bsxfun(@times,fft(val),exp(1i*shift*pi).^waveNum));

if ( realf )
    val = real(val);
end

% val(m/2,:) = mean(val(m/2,:));
% val(m,:) = mean(val(m,:));

% Create the spherefun
f = spherefun(val(m/2:m,:));

end

