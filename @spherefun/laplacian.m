function f = laplacian(f) 
%LAPLACIAN   Scalar laplacian of a SPHEREFUN.
%   L = LAPLACIAN(F).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = laplacianE(f);

end

% Compute as d2fdth2 + 1./sin(th)(cos(th) dfdth + 1./sin(th) d2fdlam2):
function f = laplacianA( f )
% We are going to work at the tech level to make things faster.
[C, D, R] = cdr( f );

% Do everything with even length columns since then no special
% modifications are required for dividing the cosine/sine series expansion.
n = length(C)+mod(length(C),2);

% Matrices for multiplying by sin/cos in coefficient space.
Msinn = .5i*spdiags(ones(n,1)*[-1,1],[-1 1],n,n);
Mcosn = .5*spdiags(ones(n,1)*[1,1],[-1 1],n,n);

% Work at the tech level to make things faster.
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

% Compute the derivatives
dCdth = diff(ctechs)/pi;
d2Cdth2 = diff(ctechs,2)/pi^2;
d2Rdlam2 = diff(rtechs,2)/pi^2;

% Calculate the C * D * R.' decomposition of 1./sin(th) d2fdlam2
C_cfs = ctechs.alias(ctechs.coeffs,n);
C1 = Msinn \ C_cfs;
R1 = d2Rdlam2.coeffs;
    
% Calculate the C * D * R.' decomposition of cos(th) dfdth
C2 = Mcosn*ctechs.alias(dCdth.coeffs,n);
R2 = rtechs.coeffs;

% Calculate the C * D * R.' decomposition of d2fdth2
C3 = ctechs.alias(d2Cdth2.coeffs,n);
R3 = rtechs.coeffs;

% Compute the g = cos(th) dfdth + 1./sin(th) d2fdlam2
g = addPieces(C1,R1,f,C2,R2,f);

% Calculate the C * D * R.' decomposition of 1./sin(th) g
if isempty(g)
    g = 0*f;
end

[C, D, R] = cdr( g );
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;
C1 = Msinn \ ctechs.alias(ctechs.coeffs,n);
R1 = rtechs.coeffs;

% Finally compute the Laplacian: d2fdth2 + g./sin(th)
f = addPieces(C1,R1,g,C3,R3,f);

f = spherefun(sample(f));
end

% Compute as 1/sin(th)^2 (sin(th)^2 d2fdth2 + sin(th)cos(th) dfdth + d2fdlam2)
function g = laplacianB( f )
% We are going to work at the tech level to make things faster.
[C, D, R] = cdr( f );

% Do everything with even length columns since then no special
% modifications are required for dividing the cosine/sine series expansion.
n = length(C)+mod(length(C),2);

% Matrices for multiplying by sin/cos in coefficient space.
Msin2n = spdiags(ones(n,1)*[-1 2 -1]/4,[-2 0 2],n,n);
Msinncosn = 0.25i*spdiags(ones(n,1)*[-1 1],[-2 2],n,n);

% Work at the tech level to make things faster.
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

% Compute the derivatives
dCdth = diff(ctechs)/pi;
d2Cdth2 = diff(ctechs,2)/pi^2;
d2Rdlam2 = diff(rtechs,2)/pi^2;

% Calculate the C * D * R.' decomposition of d2fdlam2
C1 = ctechs.alias(ctechs.coeffs,n);
R1 = d2Rdlam2.coeffs;
    
% Calculate the C * D * R.' decomposition of sin(th)cos(th) dfdth
C2 = Msinncosn*ctechs.alias(dCdth.coeffs,n);
R2 = rtechs.coeffs;

% Calculate the C * D * R.' decomposition of sin(th)^2 d2fdth2
C3 = Msin2n*ctechs.alias(d2Cdth2.coeffs,n);
R3 = rtechs.coeffs;

% Add the pieces
g = addPieces(C1,R1,f,C2,R2,f);
[C, D, R] = cdr( g );
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;
C1 = ctechs.coeffs;
R1 = rtechs.coeffs;

g = addPieces(C3,R3,f,C1,R1,g);

% Calculate the C * D * R.' decomposition of 1./sin(th)^2 g
[C, D, R] = cdr( g );
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

n = length(C)+mod(length(C),2);
Msin2n = spdiags(ones(n,1)*[-1 2 -1]/4,[-2 0 2],n,n);

C1 = Msin2n \ ctechs.alias(ctechs.coeffs,n);
R1 = rtechs.coeffs;

c1techs = real(trigtech({'',C1}));
g.cols.funs{1}.onefun = c1techs;
r1techs = real(trigtech({'',R1}));
g.rows.funs{1}.onefun = r1techs;
g.cols.pointValues = feval(c1techs,[-1;1]);
g.rows.pointValues = feval(r1techs,[-1;1]); 

end 

% Compute as d2fdth2 + cos(th)./sin(th) dfdth + 1./sin(th)^2*d2fdlam2
% THIS WAY DOES NOT WORK SINCE 1./sin(th)^2*d2fdlam2 MAY NOT BE DEFINED AT
% THE POLES.
function g = laplacianC( f )
% We are going to work at the tech level to make things faster.
[C, D, R] = cdr( f );

% Do everything with even length columns since then no special
% modifications are required for dividing the cosine/sine series expansion.
n = length(C)+mod(length(C),2);

% Matrices for multiplying by sin/cos in coefficient space.
Msinn = .5i*spdiags(ones(n,1)*[-1,1],[-1 1],n,n);
Mcosn = .5*spdiags(ones(n,1)*[1,1],[-1 1],n,n);
Msin2n = spdiags(ones(n,1)*[-1 2 -1]/4,[-2 0 2],n,n);

% Work at the tech level to make things faster.
ctechs = C.funs{1}.onefun;
rtechs = R.funs{1}.onefun;

% Compute the derivatives
dCdth = diff(ctechs)/pi;
d2Cdth2 = diff(ctechs,2)/pi^2;
d2Rdlam2 = diff(rtechs,2)/pi^2;

% Calculate the C * D * R.' decomposition of 1/sin(th).^2 d2fdlam2
C_cfs = ctechs.alias(ctechs.coeffs,n);
C1 = Msin2n \ C_cfs;
R1 = d2Rdlam2.coeffs;
    
% Calculate the C * D * R.' decomposition of cos(th)./sin(th) dfdth
C2 = Mcosn*(Msinn\ctechs.alias(dCdth.coeffs,n));
R2 = rtechs.coeffs;

% Calculate the C * D * R.' decomposition of d2fdth2
C3 = ctechs.alias(d2Cdth2.coeffs,n);

C2 = C2 + C3;

% Add the pieces
g = addPieces(C1,R1,f,C2,R2,f);

g.cols = real(g.cols);

end 

% Slow way and gives slightly higher errors than A.
function g = laplacianD( f )

% Slow way is divergence of the gradient.
[dfdx,dfdy,dfdz] = gradient(f);
[m,n] = length(f);
% Ensure n is even
n = n + mod(n,2);
% g = compose(compose(diff(dfdx,1),@plus,diff(dfdy,2)),@plus,diff(dfdz,3));
g = spherefun( sample(diff(dfdx,1),m,n/2) + ...
    sample(diff(dfdy,2),m,n/2) + sample(diff(dfdz,3),m,n/2) );

end

% Similar to D, but adds seconds derivatives together instead of sampling:
function g = laplacianE(f)
%LAPLACIANE   Do ?.

fxx = diff(f,1,2); 
fyy = diff(f,2,2); 
fzz = diff(f,3,2);

[mxx,nxx] = length(fxx);
[myy,nyy] = length(fyy);
[mzz,nzz] = length(fzz);

m = max([mxx myy mzz]);
n = max([nxx nyy nzz]);

% Ensure m and n are even:
m = m + mod(m,2);
n = n + mod(n,2);

% g = diff(f,1,2) + diff(f,2,2) + diff(f,3,2);
g = spherefun(sample(fxx, m, n/2) + sample(fyy, m, n/2) + sample(fzz, m, n/2));

end

function f = addPieces(C1,R1,f1,C2,R2,f2)

if all( R1(:,1)==0 )
    f1.pivotValues = f1.pivotValues(2:end);
    C1 = C1(:,2:end);
    R1 = R1(:,2:end);
    f1.idxPlus = f1.idxPlus(2:end)-1;
    f1.idxMinus = f1.idxMinus-1;
end
c1techs = real(trigtech({'',C1}));
f1.cols.funs{1}.onefun = c1techs;
r1techs = real(trigtech({'',R1}));
f1.rows.funs{1}.onefun = r1techs;

if all( R2(:,1)==0 )
    f2.pivotValues = f2.pivotValues(2:end);
    C2 = C2(:,2:end);
    R2 = R2(:,2:end);
    f2.idxPlus = f2.idxPlus(2:end)-1;
    f2.idxMinus = f2.idxMinus-1;
end
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

% f = spherefun(sample(f));

end
