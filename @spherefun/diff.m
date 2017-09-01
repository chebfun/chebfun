function f = diff(f, varargin)
%DIFF   Tangential derivative of a spherefun in Cartesian coordinates.
%   F = DIFF(F) computes the first tangential derivative of F with respect to x.
%   This is the projection of the surface gradient of f in the x-direction.
%
%   F = DIFF(F, DIM) computes the first tangential derivative of F. If DIM = 1,
%   the tangential derivative is taken in the x-direction. If DIM = 2, the
%   tangential derivative is taken in the y-direction and if DIM = 3, the
%   tangential derivative is taken in the z-direction.
%
%   F = DIFF(F, DIM, K) computes the Kth tangential derivatives of F in the
%   variable given by DIM.
%
% See also GRADIENT, LAPLACIAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) )
    return
end 

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

if ( dim ~= 1 && dim ~= 2 && dim ~= 3 )
    error('SPHEREFUN:DIFF:DIM', 'Unrecognized coordinate dimension.');
end

if ( abs(K - round(K)) > eps )
    error('SPHEREFUN:DIFF:DIFFORDER', 'Fractional derivatives not allowed.')
end
K = round(K);

% Implement higher derivatives as repeated (iterated) differentiation:
for j = 1:K
    f = onediff(f, dim);
end

end

function f = onediff(f, dim)
%ONEDIFF   Compute one derivative of f for the given dimension.

% [TODO]: This code will not work for complex valued spherefuns, if we ever
% allow them.

% Simplify f to avoid any extra work.
f = simplify(f);

% We are going to work at the tech level to make things faster.
[C, ~, R] = cdr(f);

% Do everything with even length columns since then no special
% modifications are required for dividing the cosine/sine series expansion.
n = length(C) + mod(length(C), 2);

% The variable coefficients in the definitions of the derivatives means
% that the length of the columns and rows will increase by one wave number
% after taking the derivatives with respect to x and y. The z derivative
% only increases the columns wave number by 1. This means we need to pad the
% coefficients with one extra zero negative and positive coefficient before
% doing the computations.
if ( dim ~= 3 )
    m = length(R) + 2;  % Pad rows
else
    m = length(R);
end
n = n + 2; % Pad columns

% Matrices for multiplying by sin/cos in coefficient space.
Msinn = .5i*spdiags(ones(n,1)*[-1,1], [-1 1], n, n);
Msinm = .5i*spdiags(ones(m,1)*[-1,1], [-1 1], m, m);
Mcosn = .5*spdiags(ones(n,1)*[1,1], [-1 1], n, n);
Mcosm = .5*spdiags(ones(m,1)*[1,1], [-1 1], m, m);

% Work at the tech level to make things faster.
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

% Alias will do the padding of the coefficients.
ctechs.coeffs = ctechs.alias(ctechs.coeffs, n);
rtechs.coeffs = rtechs.alias(rtechs.coeffs, m);

% Compute the derivatives
dCdth = diff(ctechs)/pi;
dRdlam = diff(rtechs)/pi;

% dx and dy involve two terms while dz is only one so we handle the cases
% separately:
if ( dim == 1 ) || ( dim == 2)
    
    if ( dim == 1 ) % x
        % Calculate the C * D * R.' decomposition of -sin(lam)./sin(th) dfdlam:
        C_cfs = ctechs.coeffs;
        C1 = Msinn \ C_cfs;
        R1 = -Msinm*dRdlam.coeffs;

        % Calculate the C * D * R.' decomposition of cos(lam)cos(th) dfdth:
        C2 = Mcosn*dCdth.coeffs;
        R2 = Mcosm*rtechs.coeffs;
        
    elseif ( dim == 2 ) % y
        % Calculate the C * D * R.' decomposition of cos(lam)./sin(th) dfdlam:
        C_cfs = ctechs.coeffs;
        C1 = Msinn \ C_cfs;
        R1 = Mcosm*dRdlam.coeffs;

        % Calculate the C * D * R.' decomposition of sin(lam)cos(th) dfdth:
        C2 = Mcosn*dCdth.coeffs;
        R2 = Msinm*rtechs.coeffs;
        
    end
    
    % Put pieces back together:
    f1 = f; 
    c1techs = real(trigtech({'',C1}));
    f1.cols.funs{1}.onefun = c1techs;
    r1techs = real(trigtech({'',R1}));
    f1.rows.funs{1}.onefun = r1techs;
    
    % Parity changes:
    temp = f1.idxPlus;
    f1.idxPlus = f1.idxMinus;
    f1.idxMinus = temp;

    f2 = f; 
    c2techs = real(trigtech({'',C2}));
    f2.cols.funs{1}.onefun = c2techs;
    r2techs = real(trigtech({'',R2}));
    f2.rows.funs{1}.onefun = r2techs;

    % Parity changes:
    temp = f2.idxPlus;
    f2.idxPlus = f2.idxMinus;
    f2.idxMinus = temp;

    % Compression plus may not preserve the expansion properties we want,
    % i.e. the odd/even parity properites. So we sample f1 and f2
    % separately and add their values together to construct a spherefun
    % from their samples.
    % [TODO]: Fix this so everything is done in coefficient space, 
    % using compression plus.  This will make the code more efficient.
    
    % When constructing from samples, m must be even.
    m = m + mod(m,2);
    f = spherefun(sample(f1,m,n/2+1) + sample(f2,m,n/2+1));    
    
else
    
    % Calculate the C * D * R.' decomposition of sin(th) dfdth:
    C1 = -Msinn*dCdth.coeffs;
    R1 = rtechs.coeffs;

    % Put pieces back together:
    c1techs = real(trigtech({'',C1}));
    f.cols.funs{1}.onefun = c1techs;
    r1techs = real(trigtech({'',R1}));
    f.rows.funs{1}.onefun = r1techs;    

    % Weird feval behavior in CHEBFUN requires this:
    f.cols.pointValues = feval(c1techs, [-1;1]);
    f.rows.pointValues = feval(r1techs, [-1;1]); 
    
end    

end
