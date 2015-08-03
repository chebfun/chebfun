function f = rotate( f, alpha, beta)
%ROTATE  Rotates a SPHEREFUN by angles (alpha,beta)
%   Y = ROTATE( F, ALPHA, BETA ) rotates f by an angle alpha in the (x,z)
%   plane and beta in the (x,y) plane.

if nargin == 1
    return;
elseif nargin == 2
    beta = 0;
end

% Determine the transformation to use
if iscolat( f )
    fr = @(lam,th) rotateColat(f,lam,th,alpha,beta);
else
    fr = @(lam,th) rotateLat(f,lam,th,alpha,beta);
end
f = spherefun(fr, f.domain);
end

function y = rotateColat(f,lam,th,a,b)
    cosa = cos(a); cosb = cos(b);
    sina = sin(a); sinb = sin(b);
    x = cos(lam).*sin(th);
    y = sin(lam).*sin(th);
    z = cos(th);
    
    u = cosb*cosa*x - sinb*y - cosb*sina*z;
    v = sinb*cosa*x + cosb*y - sinb*sina*z;
    w = sina*x + cosa*z;
    
    [lam,th] = cart2sph(u,v,w);
    
    y = feval(f,lam,pi/2-th);
end
    
function y = rotateLat(f,lam,th,a,b)
    cosa = cos(a); cosb = cos(b);
    sina = sin(a); sinb = sin(b);
    x = cos(lam).*cos(th);
    y = sin(lam).*cos(th);
    z = sin(th);
    
    u = cosb*cosa*x - sinb*y - cosb*sina*z;
    v = sinb*cosa*x + cosb*y - sinb*sina*z;
    w = sina*x + cosa*z;
    
    [lam,th] = cart2sph(u,v,w);
    
    y = feval(f,lam,th);
end
