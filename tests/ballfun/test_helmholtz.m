function pass = test_helmholtz( )
% Test the Helmholtz solver with Dirichlet BC
S = [39,40,41];

%% POISSON WITH DIRICHLET BOUNDARY CONDITIONS

% Example 1:
f = cheb.galleryballfun('zero',S);
exact = @(r, lam, th) 1;
bc = @(lam, th) exact(1, lam, th);  % Dirichlet condition
u = helmholtz(f, 0, bc);
exact = ballfun(exact,S);
pass(1) = isequal(u,exact);

% Example 2:
f = ballfun(@(r, lam, th) 6,S);
exact = @(r, lam, th) r.^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc);
exact = ballfun(exact,S);
pass(2) = isequal(u,exact);

% Example 3:
f = ballfun(@(r, lam, th) 4,S);
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc);
exact = ballfun(exact,S);
pass(3) = isequal(u,exact);

% Example 4:
f = cheb.galleryballfun('zero',S);       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = helmholtz(f, 0, bc);
exact = ballfun(exact,S);
pass(4) = isequal(u,exact);
    
% Example 5:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)),S);
exact = @(r, lam, th) cos( r.^2.*sin(th).^2.*cos(lam).^2 );
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, 0, bc);
exact = ballfun(exact,S);
pass(5) = isequal(u,exact);

S = [50,50,50];
zero = cheb.galleryballfun("zero",S);

% Example 6:
g = @(r,lam,th)r.^2.*sin(th).^2.*exp(2*1i*lam);
f = ballfun(g,S);
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero,0,bc);
pass(6) = isequal(laplace(f),zero) && isequal(u,f);

% Example 7:
g = @(r,lam,th)r.^5.*exp(3*1i*lam).*sin(th).^3.*(9*cos(th).^2-1);
f = ballfun(g,S);
bc = @(lam,th)g(1,lam,th);
u = helmholtz(zero,0,bc);
pass(7) = isequal(laplace(f),zero) && isequal(u,f);

%% HELMHOLTZ WITH K = 2

K = 2;

% Example 8:
f = ballfun(@(r,lam,th)4+0*r,S);
exact = @(r, lam, th)1+0*r;
bc = @(lam, th) exact(1, lam, th);  % Dirichlet condition
u = helmholtz(f, K, bc);
exact = ballfun(exact,S);
pass(8) = isequal(u,exact);

% Example 9:
f = ballfun(@(r, lam, th) 6 + 4*r.^2,S);
exact = @(r, lam, th) r.^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc);
exact = ballfun(exact,S);
pass(9) = isequal(u,exact);

% Example 10:
f = ballfun(@(r, lam, th) 4 + 4*r.^2.*sin(th).^2,S);
exact = @(r, lam, th) r.^2.*sin(th).^2;
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc);
exact = ballfun(exact,S);
pass(10) = isequal(u,exact);

% Example 11:
f = ballfun(@(r,lam,th)4*r.*sin(lam).*sin(th),S);       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = helmholtz(f, K, bc);
exact = ballfun(exact,S);
pass(11) = isequal(u,exact);
    
% Example 12:
f = ballfun(@(r, lam, th) -2*(2*r.^2.*cos(lam).^2.*cos(r.^2.*cos(lam).^2.*sin(th).^2).*sin(th).^2 + ...
    sin(r.^2.*cos(lam).^2.*sin(th).^2)) + 4*cos(r.^2.*sin(th).^2.*cos(lam).^2),S);
exact = @(r, lam, th) cos(r.^2.*sin(th).^2.*cos(lam).^2);
bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
u = helmholtz(f, K, bc);
exact = ballfun(exact,S);
pass(12) = isequal(u,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end
