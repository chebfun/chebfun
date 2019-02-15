function pass = test_helmholtz( ) 
% Test the Helmholtz solver with Dirichlet BC

%% POISSON WITH DIRICHLET BOUNDARY CONDITIONS

% Example 1:
f = ballfun(@(x,y,z)0);
exact = @(r, lam, th) 1;
bc = @(lam, th) exact(1, lam, th);  % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(1) = isequal(u,exact);

% Example 2:
f = ballfun(@(r, lam, th) 6);
exact = @(r, lam, th) r.^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(2) = isequal(u,exact);

% Example 3:
f = ballfun(@(r, lam, th) 4);
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(3) = isequal(u,exact);

% Example 4:
f = ballfun(@(x,y,z)0);       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(4) = isequal(u,exact);

% Example 5:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)),'spherical');
exact = @(r, lam, th) cos( r.^2.*sin(th).^2.*cos(lam).^2 );
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact,'spherical');
pass(5) = isequal(u,exact);

zero = ballfun(@(x,y,z)0);

% Example 6:
g = @(r,lam,th)r.^2.*sin(th).^2.*exp(2*1i*lam);
f = ballfun(g,'spherical');
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero, 0, bc, 50, 50, 50);
pass(6) = isequal(laplacian(u),zero) && isequal(u,f);

% Example 7:
g = @(r,lam,th)r.^5.*exp(3*1i*lam).*sin(th).^3.*(9*cos(th).^2-1);
f = ballfun(g,'spherical');
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero, 0, bc, 50, 50, 50);
pass(7) = isequal(laplacian(u),zero) && isequal(u,f);

%% HELMHOLTZ WITH K = 2 AND DIRICHLET BC

K = 2;

% Example 8:
f = ballfun(@(r,lam,th)4+0*r,'spherical');
exact = @(r, lam, th)1+0*r;
bc = @(lam, th) exact(1, lam, th);  % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(8) = isequal(u,exact);

% Example 9:
f = ballfun(@(r, lam, th) 6 + 4*r.^2,'spherical');
exact = @(r, lam, th) r.^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(9) = isequal(u,exact);

% Example 10:
f = ballfun(@(r, lam, th) 4 + 4*r.^2.*sin(th).^2,'spherical');
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(10) = isequal(u,exact);

% Example 11:
f = ballfun(@(r,lam,th)4*r.*sin(lam).*sin(th),'spherical');       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(11) = isequal(u,exact);
    
% Example 12:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2),'spherical');
exact = @(r, lam, th) cos(r.^2.*sin(th).^2.*cos(lam).^2);
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact,'spherical');
pass(12) = isequal(u,exact);

%% POISSON WITH NEUMANN BC

% Example 13:
exact = ballfun(@(r, lam, th)1,'spherical');
f = laplacian(exact);
bc = @(lam, th) 0;  % Neumann condition
u = helmholtz(f,0,bc,39,40,41,'neumann');
pass(13) = isequal(diff(u,1,'spherical'),diff(exact,1,'spherical')) ...
        && isequal(diff(u,2,'spherical'),diff(exact,2,'spherical')) ...
        && isequal(diff(u,3,'spherical'),diff(exact,3,'spherical'));

% Example 14:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'spherical');
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = helmholtz(f,0,bc,39,40,41,'neumann');
pass(14) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 15:
exact = ballfun(@(r, lam, th)r.^3.*sin(th).^3.*cos(lam),'spherical');
f = laplacian(exact);
bc = @(lam, th) 3*sin(th).^3.*cos(lam);  % Neumann condition
u = helmholtz(f,0,bc,39,40,41,'neumann');
pass(15) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 16:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2.*cos(lam).^2,'spherical');
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2.*cos(lam).^2;  % Neumann condition
u = helmholtz(f,0,bc,39,40,41,'neumann');
pass(16) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 17:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,40,54,41,'neumann');
pass(17) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 18:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,41,54,41,'neumann');
pass(18) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 19:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,42,54,41,'neumann');
pass(19) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 20:
exact = ballfun(@(r,lam,th)(r.*sin(th).*cos(lam).*r.*cos(th)).^2,'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,54,41,'neumann');
pass(20) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 21:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,54,41,'neumann');
pass(21) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 22:
exact = ballfun(@(r,lam,th)cos((r.*sin(th).*cos(lam)).^3),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,54,41,'neumann');
pass(22) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 23:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,43,54,42,'neumann');
pass(23) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 24:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,43,54,43,'neumann');
pass(24) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 25:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,43,54,44,'neumann');
pass(25) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 26:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz(f,0,bc,42,42,42,'neumann');
pass(26) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 27:
exact = ballfun(@(r,lam,th)(r.*sin(th).*cos(lam).*r.*cos(th)).^2,'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,42,42,'neumann');
pass(27) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 28:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,42,42,'neumann');
pass(28) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 29:
exact = ballfun(@(r,lam,th)cos((r.*sin(th).*cos(lam)).^3),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,42,42,'neumann');
pass(29) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

%% HELMHOLTZ WITH K = 2  AND NEUMANN BC

K = 2;

% Example 30:
exact = ballfun(@(r, lam, th)1,'spherical');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 0;  % Neumann condition
u = helmholtz(f, K, bc1, 38, 17, 24,'neumann');
pass(30) = isequal(u,exact);

% Example 31:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'spherical');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = helmholtz(f, K, bc1,38,17,24,'neumann');
pass(31) = isequal(u,exact);

% Example 32:
exact = ballfun(@(r, lam, th)r.^4.*sin(th).^2,'spherical');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 4*sin(th).^2;  % Neumann condition
u = helmholtz(f, K, bc1,38,17,24,'neumann');
pass(32) = isequal(u,exact);

%% M = N = P

% Example 33:
g = @(r,lam,th)r.^5.*exp(3*1i*lam).*sin(th).^3.*(9*cos(th).^2-1);
f = ballfun(g,'spherical');
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero, 0, bc, 50);
pass(33) = isequal(laplacian(u),zero) && isequal(u,f);

% Example 34:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),'spherical');
f = laplacian(exact);
bc = diff(exact,1,'spherical');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz(f,0,bc,42,'neumann');
pass(34) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

if (nargout > 0)
    pass = all(pass(:));
end
end
