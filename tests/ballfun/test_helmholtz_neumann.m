function pass = test_helmholtz_neumann( )

eps = 1e-10;

%% NEUMANN BOUNDARY CONDITIONS
S = [39,40,41];
n = S(2); p = S(3);

% Example 1:
f = ballfun(@(r, lam, th)cos(lam),S);
bc = @(lam, th) sin(lam).*cos(th);   % Neumann condition
u = helmholtz_neumann(f,0,bc);
U = u.coeffs;
pass(1) = abs(U(1,floor(n/2)+1,floor(p/2)+1))<eps;

% Example 2:
f = ballfun(@(r, lam, th)r.*cos(lam).*sin(th),S);
bc = @(lam, th) 3*sin(lam).*cos(th);   % Neumann condition
u = helmholtz_neumann(f,0,bc);
U = u.coeffs;
pass(2) = abs(U(1,floor(n/2)+1,floor(p/2)+1))<eps;

% Example 3:
f = ballfun(@(r, lam, th)r.^2.*sin(th).^2,S);
bc = @(lam, th) 2*sin(th).^2;   % Neumann condition
u = helmholtz_neumann(f,0,bc);
U = u.coeffs;
pass(3) = abs(U(1,floor(n/2)+1,floor(p/2)+1))<eps;

% Example 4:
exact = ballfun(@(r, lam, th)1,S);
f = laplacian(exact);
bc = @(lam, th) 0;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(4) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 5:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,S);
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(5) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 6:
exact = @(r, lam, th) r.^2.*sin(th).^2.*cos(lam);
f = ballfun(@(r, lam, th) 3*cos(lam),S);
bc = @(lam, th) 2*sin(th).^2.*cos(lam);  % Neumann condition
u = helmholtz_neumann(f,0,bc);
exact = ballfun(exact,S);
pass(6) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 7:
exact = ballfun(@(r, lam, th)r.^3.*sin(th).^5,S);
f = laplacian(exact);
bc = @(lam, th) 3*sin(th).^5;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(7) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 8:
exact = ballfun(@(r, lam, th)r.^3.*sin(th).^3.*cos(lam),S);
f = laplacian(exact);
bc = @(lam, th) 3*sin(th).^3.*cos(lam);  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(8) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 9:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2.*cos(lam).^2,S);
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2.*cos(lam).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(9) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 10:
exact = ballfun(@(r, lam, th)r.^3.*sin(th).^3.*cos(lam).^2,S);
f = laplacian(exact);
bc = @(lam, th) 3*sin(th).^3.*cos(lam).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(10) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

S = [40,54,41];

% Example 11:
exact = ballfun(@(x,y,z)y.^2,'cart',S);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(11) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

S = [41,54,41];

% Example 12:
exact = ballfun(@(x,y,z)y.^2,'cart',S);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(12) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
S = [42,54,41];

% Example 13:
exact = ballfun(@(x,y,z)y.^2,'cart',S);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(13) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 14:
exact = ballfun(@(r,lam,th)(r.*sin(th).*cos(lam).*r.*cos(th)).^2,S);
f = laplacian(exact);
bc = diff(exact,1);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc);
pass(14) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 15:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),S);
f = laplacian(exact);
bc = diff(exact,1);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc);
pass(15) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 16:
exact = ballfun(@(r,lam,th)cos((r.*sin(th).*cos(lam)).^3),S);
f = laplacian(exact);
bc = diff(exact,1);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc);
pass(16) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

S = [43,54,42];

% Example 17:
exact = ballfun(@(x,y,z)y.^2,'cart',S);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(17) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

S = [43,54,43];

% Example 18:
exact = ballfun(@(x,y,z)y.^2,'cart',S);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(18) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

S = [43,54,44];

% Example 19:
exact = ballfun(@(x,y,z)y.^2,'cart',S);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(19) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
S = [42,42,42];

% Example 20:
exact = ballfun(@(x,y,z)y.^2,'cart',S);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc);
pass(20) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 21:
exact = ballfun(@(r,lam,th)(r.*sin(th).*cos(lam).*r.*cos(th)).^2,S);
f = laplacian(exact);
bc = diff(exact,1);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc);
pass(21) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 22:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),S);
f = laplacian(exact);
bc = diff(exact,1);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc);
pass(22) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 23:
exact = ballfun(@(r,lam,th)cos((r.*sin(th).*cos(lam)).^3),S);
f = laplacian(exact);
bc = diff(exact,1);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc);
pass(23) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

%% HELMHOLTZ WITH K = 2  AND NEUMANN BC

S = [38,17,24];
K = 2;

% Example 24:
exact = ballfun(@(r, lam, th)1,S);
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 0;  % Neumann condition
u = helmholtz_neumann(f, K, bc1);
pass(24) = isequal(u,exact);

% Example 25:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,S);
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = helmholtz_neumann(f, K, bc1);
pass(25) = isequal(u,exact);

% Example 26:
exact = ballfun(@(r, lam, th) r.^2.*sin(th).^2.*cos(lam),S);
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 2*sin(th).^2.*cos(lam);  % Neumann condition
u = helmholtz_neumann(f, K, bc1);
pass(26) = isequal(u,exact);

% Example 27:
exact = ballfun(@(r, lam, th)r.^4.*sin(th).^2,S);
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 4*sin(th).^2;  % Neumann condition
u = helmholtz_neumann(f, K, bc1);
pass(27) = isequal(u,exact);

% Example 28:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^4.*cos(lam),S);
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 2*sin(th).^4.*cos(lam);  % Neumann condition
u = helmholtz_neumann(f, K, bc1);
pass(28) = isequal(u,exact);

% Example 29:
exact = ballfun(@(r, lam, th)r.^3.*sin(th).^5.*cos(lam).^2,S);
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 3*sin(th).^5.*cos(lam).^2;  % Neumann condition
u = helmholtz_neumann(f, K, bc1);
pass(29) = isequal(u,exact);

% Example 30:
exact = ballfun(@(r, lam, th)r.^3.*sin(th).^3.*cos(lam).^2,S);
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 3*sin(th).^3.*cos(lam).^2;  % Neumann condition
u = helmholtz_neumann(f, K, bc1);
pass(30) = isequal(u,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end
