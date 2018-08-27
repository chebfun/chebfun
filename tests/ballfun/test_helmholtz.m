function pass = test_helmholtz( ) 
% Test the Helmholtz solver with Dirichlet BC

%% POISSON WITH DIRICHLET BOUNDARY CONDITIONS

% Example 1:
f = cheb.galleryballfun('zero');
exact = @(r, lam, th) 1;
bc = @(lam, th) exact(1, lam, th);  % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact);
pass(1) = isequal(u,exact);

% Example 2:
f = ballfun(@(r, lam, th) 6);
exact = @(r, lam, th) r.^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact);
pass(2) = isequal(u,exact);

% Example 3:
f = ballfun(@(r, lam, th) 4);
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact);
pass(3) = isequal(u,exact);

% Example 4:
f = cheb.galleryballfun('zero');       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact);
pass(4) = isequal(u,exact);

% Example 5:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)));
exact = @(r, lam, th) cos( r.^2.*sin(th).^2.*cos(lam).^2 );
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc, 39, 40, 41);
exact = ballfun(exact);
pass(5) = isequal(u,exact);

zero = cheb.galleryballfun("zero");

% Example 6:
g = @(r,lam,th)r.^2.*sin(th).^2.*exp(2*1i*lam);
f = ballfun(g);
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero, 0, bc, 50, 50, 50);
pass(6) = isequal(laplacian(u),zero) && isequal(u,f);

% Example 7:
g = @(r,lam,th)r.^5.*exp(3*1i*lam).*sin(th).^3.*(9*cos(th).^2-1);
f = ballfun(g);
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero, 0, bc, 50, 50, 50);
pass(7) = isequal(laplacian(u),zero) && isequal(u,f);

%% HELMHOLTZ WITH K = 2

K = 2;

% Example 8:
f = ballfun(@(r,lam,th)4+0*r);
exact = @(r, lam, th)1+0*r;
bc = @(lam, th) exact(1, lam, th);  % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact);
pass(8) = isequal(u,exact);

% Example 9:
f = ballfun(@(r, lam, th) 6 + 4*r.^2);
exact = @(r, lam, th) r.^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact);
pass(9) = isequal(u,exact);

% Example 10:
f = ballfun(@(r, lam, th) 4 + 4*r.^2.*sin(th).^2);
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact);
pass(10) = isequal(u,exact);

% Example 11:
f = ballfun(@(r,lam,th)4*r.*sin(lam).*sin(th));       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact);
pass(11) = isequal(u,exact);
    
% Example 12:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2));
exact = @(r, lam, th) cos(r.^2.*sin(th).^2.*cos(lam).^2);
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc, 50, 50, 50);
exact = ballfun(exact);
pass(12) = isequal(u,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end
