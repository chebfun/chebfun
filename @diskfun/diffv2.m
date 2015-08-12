function f = diff( f, varargin )
% DIFF    Derivative of a diskfun in Cartesian coordinates.
%
%  F = DIFF( F ) computes the first derivative of F with respect to x.
%
%  F = DIFF( F, DIM )  computes the first derivative of F. If DIM = 1, the
%  derivative is taken in the x-direction. If DIM = 2, the derivative
%  is taken in the y-direction.
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
    error('DISKFUN:DIFF:ORDER', 'Only first derivatives currently allowed.');
end

if ( dim ~= 1 && dim ~= 2 )
    error('DISKFUN:DIFF:DIM', 'Unrecognized coordinate dimension');
end

if ( abs( K - round(K) ) > eps )
    error('DISKFUN:DIFF:DIFFORDER', 'Fractional derivatives not allowed')
end
K = round( K );

% Track whether f is real.
realf = isreal(f);

% We are going to work at the tech level to make things faster.
[cols, D, rows] = cdr( f );

coltechs = cols.funs{1}.onefun; %r [-1, 1]
rowtechs = rows.funs{1}.onefun; %theta [-pi,pi]

% Compute the derivatives
d_coltechs = diff(coltechs);
d_rowtechs = diff(rowtechs)/pi; %COV 

% We will do everyting in value space on [-1,1] at the tech level.
%note: chebtech2 doesn't store values?

% Evaluate at even chebyshev grid points in r so pole is not evaluated
m = length(cols);
n = length(rows);

%alias the coeffs over an even grid (m+1) and then use coeffs2val

coltechs2= chebtech2.alias(coltechs.coeffs, m+1);
d_coltechs2=chebtech2.alias(d_coltechs.coeffs, m+1); 

colv=chebtech2.coeffs2vals(coltechs2);
dcolv=chebtech2.coeffs2vals(d_coltechs2); 

% r and theta at appropriate grids
th = trigpts(n,[-1,1]).' ;
r = chebpts(m+1,[-1,1]);

% Evalute rows and the derivatives at the grid
rowv = rowtechs.values.';
drowv = d_rowtechs.values.';

if ( dim == 1 )            % x
  
     val = repmat(cos(th),m+1,1).*(dcolv*D*rowv) -...
         (1./r)*sin(th).*(colv*D*drowv);
        
        
else ( dim == 2 )         % y
    
     val = repmat(sin(th),m+1,1).*(dcolv*D*rowv) +...
         (1./r)*cos(th).*(colv*D*drowv);

end

%
% Shift back to regular grid points in r and then build a diskfun




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

