function pass = test_helmholtz( pref ) 
% Test the Helmholtz solver with Dirichlet BC

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e7*pref.techPrefs.chebfuneps;

%% POISSON WITH DIRICHLET BOUNDARY CONDITIONS

% Example 1:
f = ballfun(@(x,y,z)0);
exact = @(r, lam, th) 1;
bc = @(lam, th) exact(1, lam, th);  % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(1) = norm(u - exact) < tol;

% Example 2:
f = ballfun(@(r, lam, th) 6);
exact = @(r, lam, th) r.^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(2) = norm(u - exact) < tol;

% Example 3:
f = ballfun(@(r, lam, th) 4);
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(3) = norm(u - exact) < tol;

% Example 4:
f = ballfun(@(x,y,z)0);       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(4) = norm(u - exact) < tol;

% Example 5:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)),'spherical');
exact = @(r, lam, th) cos( r.^2.*sin(th).^2.*cos(lam).^2 );
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(5) = norm(u - exact) < tol;

zero = ballfun(0);

% Example 6:
g = @(r,lam,th)r.^2.*sin(th).^2.*exp(2*1i*lam);
f = ballfun(g,'spherical');
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero, 0, bc, 50, 50, 50);
pass(6) = norm(laplacian(u)) < tol && norm(u - f) < tol;

% Example 7:
g = @(r,lam,th)r.^5.*exp(3*1i*lam).*sin(th).^3.*(9*cos(th).^2-1);
f = ballfun(g,'spherical');
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero, 0, bc, 50, 50, 50);
pass(7) = norm(laplacian(u)) < tol && norm(u - f) < tol;

%% HELMHOLTZ WITH K = 2 AND DIRICHLET BC

K = 2;

% Example 8:
f = ballfun(@(r,lam,th)4+0*r,'spherical');
exact = @(r, lam, th)1+0*r;
bc = @(lam, th) exact(1, lam, th);  % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(8) = norm(u - exact) < tol;

% Example 9:
f = ballfun(@(r, lam, th) 6 + 4*r.^2,'spherical');
exact = @(r, lam, th) r.^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(9) = norm(u - exact) < tol;

% Example 10:
f = ballfun(@(r, lam, th) 4 + 4*r.^2.*sin(th).^2,'spherical');
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(10) = norm(u - exact) < tol;

% Example 11:
f = ballfun(@(r,lam,th)4*r.*sin(lam).*sin(th),'spherical');       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(11) = norm(u - exact) < tol;
    
% Example 12:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2),'spherical');
exact = @(r, lam, th) cos(r.^2.*sin(th).^2.*cos(lam).^2);
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(12) = norm(u - exact) < tol;

%% POISSON WITH NEUMANN BC

% Example 13:
exact = ballfun(@(r, lam, th)1,'spherical');
f = laplacian(exact);
bc = @(lam, th) 0;  % Neumann condition
u = helmholtz(f,0,bc,39,40,41,'neumann');
pass(13) = norm(diff(u,1,'spherical') - diff(exact,1,'spherical')) < tol ...
        && norm(diff(u,2,'spherical') - diff(exact,2,'spherical')) < tol ...
        && norm(diff(u,3,'spherical') - diff(exact,3,'spherical')) < tol;

% Example 14:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'spherical');
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = helmholtz(f,0,bc,39,40,41,'neumann');
pass(14) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 15:
exact = ballfun(@(r, lam, th)r.^3.*sin(th).^3.*cos(lam),'spherical');
f = laplacian(exact);
bc = @(lam, th) 3*sin(th).^3.*cos(lam);  % Neumann condition
u = helmholtz(f,0,bc,39,40,41,'neumann');
pass(15) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;
    
% Example 16:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2.*cos(lam).^2,'spherical');
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2.*cos(lam).^2;  % Neumann condition
u = helmholtz(f,0,bc,39,40,41,'neumann');
pass(16) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 17:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,40,54,41,'neumann');
pass(17) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 18:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,41,54,41,'neumann');
pass(18) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;
    
% Example 19:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,42,54,41,'neumann');
pass(19) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;
    
% Example 20:
exact = ballfun(@(r,lam,th)(r.*sin(th).*cos(lam).*r.*cos(th)).^2,'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,54,41,'neumann');
pass(20) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 21:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,54,41,'neumann');
pass(21) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 22:
exact = ballfun(@(r,lam,th)cos((r.*sin(th).*cos(lam)).^3),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,54,41,'neumann');
pass(22) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 23:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,43,54,42,'neumann');
pass(23) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 24:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,43,54,43,'neumann');
pass(24) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 25:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,43,54,44,'neumann');
pass(25) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;
    
% Example 26:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,42,42,42,'neumann');
pass(26) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 27:
exact = ballfun(@(r,lam,th)(r.*sin(th).*cos(lam).*r.*cos(th)).^2,'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,42,42,'neumann');
pass(27) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 28:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,42,42,'neumann');
pass(28) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 29:
exact = ballfun(@(r,lam,th)cos((r.*sin(th).*cos(lam)).^3),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,42,42,'neumann');
pass(29) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

%% HELMHOLTZ WITH K = 2  AND NEUMANN BC

K = 2;

% Example 30:
exact = ballfun(@(r, lam, th)1,'spherical');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 0;  % Neumann condition
u = helmholtz(f, K, bc1, 38, 17, 24,'neumann');
pass(30) = norm(u - exact) < tol;

% Example 31:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'spherical');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = helmholtz(f, K, bc1,38,17,24,'neumann');
pass(31) = norm(u - exact) < tol;

% Example 32:
exact = ballfun(@(r, lam, th)r.^4.*sin(th).^2,'spherical');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 4*sin(th).^2;  % Neumann condition
u = helmholtz(f, K, bc1,38,17,24,'neumann');
pass(32) = norm(u - exact) < tol;

%% M = N = P

% Example 33:
g = @(r,lam,th)r.^5.*exp(3*1i*lam).*sin(th).^3.*(9*cos(th).^2-1);
f = ballfun(g,'spherical');
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero, 0, bc, 50);
pass(33) = norm(laplacian(u)) < tol && norm(u - f) < tol;

% Example 34:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,'neumann');
pass(34) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

%% Spherefun boundary conditions

% Example 35:
exact = ballfun(@(x,y,z)(x.^2+y.^2+z.^2).*(x.^2+y.^2));
f = laplacian(exact) + 4*exact;
bc = spherefun(@(x,y,z)4*(x.^2+y.^2));  % Neumann condition
u = helmholtz(f, 2, bc, 40, 'neumann');
pass(35) = norm(u - exact) < tol;

% Example 36:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = spherefun(@(lam, th) 2*(sin(th).*sin(lam)).^2);  % Neumann condition
u = helmholtz(f,0,bc,43,54,44,'neumann');
pass(36) = norm(diff(u,1) - diff(exact,1)) < tol ...
        && norm(diff(u,2) - diff(exact,2)) < tol ...
        && norm(diff(u,3) - diff(exact,3)) < tol;

% Example 37:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)),'spherical');
exact = @(r, lam, th) cos( r.^2.*sin(th).^2.*cos(lam).^2 );
bc = spherefun(@(lam,th)cos(sin(th).^2.*cos(lam).^2 ));   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(37) = norm(u - exact) < tol;

% Example 38:
f = ballfun(@(r, lam, th) 4);
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = spherefun(@(lam, th) sin(th).^2);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(38) = norm(u - exact) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
