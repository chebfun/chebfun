function pass = test_poisson()

%% Poisson with Dirichlet BC
% Example 1:
f = cheb.galleryballfun('zero');       % Forcing function
exact = @(r, lam, th) r.*sin(lam).*sin(th);
bc = @(lam, th) exact(1, th, lam);           % Dirichlet condition around r=1
u = poisson(f, bc, 39, 40, 41);
exact = ballfun(exact,'polar');
pass(1) = isequal(u,exact);

%% Poisson with Neumann BC
% Example 2
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'polar');
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = poisson(f, bc, 39, 40, 41, 'neumann');
pass(2) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

if (nargout > 0)
    pass = all(pass(:));
end
end

