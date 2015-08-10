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

% TODO: This code will not work for complex valued spherefuns, if we ever
% allow them.
% realf = isreal(f);

% We are going to work at the tech level to make things faster.
[C, D, R] = cdr( f );

% Do everything with even length columns since then no special
% modifications are required for dividing the cosine/sine series expansion.
n = length(C)+mod(length(C),2);
m = length(R);

% Matrices for multiplying by sin/cos in coefficient space.
Msinn = .5i*spdiags(ones(n,1)*[-1,1],[-1 1],n,n);
Msinm = .5i*spdiags(ones(m,1)*[-1,1],[-1 1],m,m);
Mcosn = .5*spdiags(ones(n,1)*[1,1],[-1 1],n,n);
Mcosm = .5*spdiags(ones(m,1)*[1,1],[-1 1],m,m);

% Work at the tech level to make things faster.
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

% Compute the derivatives
dCdth = diff(ctechs)/pi;
dRdlam = diff(rtechs)/pi;

if ( dim == 1 )            % x
    % Calculate the C * D * R.' decomposition of -sin(lam)./sin(th) dfdlam
    C_cfs = ctechs.alias(ctechs.coeffs,n);
    C1 = Msinn \ C_cfs;
    R1 = -Msinm*dRdlam.coeffs;
    
    % Calculate the C * D * R.' decomposition of cos(lam)cos(th) dfdth
    C2 = Mcosn*ctechs.alias(dCdth.coeffs,n);
    R2 = Mcosm*rtechs.coeffs;
elseif ( dim == 2 )         % y
    % Calculate the C * D * R.' decomposition of cos(lam)./sin(th) dfdlam
    C_cfs = ctechs.alias(ctechs.coeffs,n);
    C1 = Msinn \ C_cfs;
    R1 = Mcosm*dRdlam.coeffs;
    
    % Calculate the C * D * R.' decomposition of sin(lam)cos(th) dfdth
    C2 = Mcosn*ctechs.alias(dCdth.coeffs,n);
    R2 = Msinm*rtechs.coeffs;
else
    % No dfdlam term;
    C1 = 0;
    R1 = 0;
    
    % Calculate the C * D * R.' decomposition of sin(th) dfdth
    C2 = -Msinn*ctechs.alias(dCdth.coeffs,n);
    R2 = rtechs.coeffs;
end
% Put pieces back together
f1 = f; 
c1techs = real(trigtech({'',C1}));
f1.cols.funs{1}.onefun = c1techs;
r1techs = real(trigtech({'',R1}));
f1.rows.funs{1}.onefun = r1techs;

f2 = f; 
c2techs = real(trigtech({'',C2}));
f2.cols.funs{1}.onefun = c2techs;
r2techs = real(trigtech({'',R2}));
f2.rows.funs{1}.onefun = r2techs;

% Weird feval behavior in chebfun requires this
f1.cols.pointValues = feval(c1techs,[-1;1]);
f1.rows.pointValues = feval(r1techs,[-1;1]); 
f2.cols.pointValues = feval(c2techs,[-1;1]);
f2.rows.pointValues = feval(r2techs,[-1;1]);

f = f1 + f2;

end

