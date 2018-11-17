function pass = test_helmholtz_neumann( pref ) 
% Test the Helmholtz solver with neumann bc. 

%% NEUMANN BOUNDARY CONDITIONS
n = 40; p = 41;

% Example 1:
f = ballfun(@(r, lam, th)r.*cos(lam).*sin(th),'polar');
bc = @(lam, th) 3*sin(lam).*cos(th);   % Neumann condition
u = helmholtz_neumann(f,0,bc,39,40,41);
U = u.coeffs;
pass(1) = abs(U(1,floor(n/2)+1,floor(p/2)+1))<eps;

% Example 2:
f = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'polar');
bc = @(lam, th) 2*sin(th).^2;   % Neumann condition
u = helmholtz_neumann(f,0,bc,39,40,41);
U = u.coeffs;
pass(2) = abs(U(1,floor(n/2)+1,floor(p/2)+1))<eps;

% Example 3:
exact = ballfun(@(r, lam, th)1,'polar');
f = laplacian(exact);
bc = @(lam, th) 0;  % Neumann condition
u = helmholtz_neumann(f,0,bc,39,40,41);
pass(3) = isequal(diff(u,1,'polar'),diff(exact,1,'polar')) ...
        && isequal(diff(u,2,'polar'),diff(exact,2,'polar')) ...
        && isequal(diff(u,3,'polar'),diff(exact,3,'polar'));

% Example 4:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'polar');
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,39,40,41);
pass(4) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 5:
exact = ballfun(@(r, lam, th)r.^3.*sin(th).^3.*cos(lam),'polar');
f = laplacian(exact);
bc = @(lam, th) 3*sin(th).^3.*cos(lam);  % Neumann condition
u = helmholtz_neumann(f,0,bc,39,40,41);
pass(5) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 6:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2.*cos(lam).^2,'polar');
f = laplacian(exact);
bc = @(lam, th) 2*sin(th).^2.*cos(lam).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,39,40,41);
pass(6) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 7:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,40,54,41);
pass(7) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 8:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,41,54,41);
pass(8) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 9:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,42,54,41);
pass(9) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 10:
exact = ballfun(@(r,lam,th)(r.*sin(th).*cos(lam).*r.*cos(th)).^2,'polar');
f = laplacian(exact);
bc = diff(exact,1,'polar');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc,42,54,41);
pass(10) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 11:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),'polar');
f = laplacian(exact);
bc = diff(exact,1,'polar');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc,42,54,41);
pass(11) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 12:
exact = ballfun(@(r,lam,th)cos((r.*sin(th).*cos(lam)).^3),'polar');
f = laplacian(exact);
bc = diff(exact,1,'polar');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc,42,54,41);
pass(12) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 13:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,43,54,42);
pass(13) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 14:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,43,54,43);
pass(14) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 15:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,43,54,44);
pass(15) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));
    
% Example 16:
exact = ballfun(@(x,y,z)y.^2);
f = laplacian(exact);
bc = @(lam, th) 2*(sin(th).*sin(lam)).^2;  % Neumann condition
u = helmholtz_neumann(f,0,bc,42,42,42);
pass(16) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 17:
exact = ballfun(@(r,lam,th)(r.*sin(th).*cos(lam).*r.*cos(th)).^2,'polar');
f = laplacian(exact);
bc = diff(exact,1,'polar');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc,42,42,42);
pass(17) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 18:
exact = ballfun(@(r,lam,th)sin((r.*sin(th).*sin(lam)).^2),'polar');
f = laplacian(exact);
bc = diff(exact,1,'polar');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc,42,42,42);
pass(18) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

% Example 19:
exact = ballfun(@(r,lam,th)cos((r.*sin(th).*cos(lam)).^3),'polar');
f = laplacian(exact);
bc = diff(exact,1,'polar');
S = size(bc);
bc = reshape(sum(bc.coeffs,1),S(2),S(3));
u = helmholtz_neumann(f,0,bc,42,42,42);
pass(19) = isequal(diff(u,1),diff(exact,1)) ...
        && isequal(diff(u,2),diff(exact,2)) ...
        && isequal(diff(u,3),diff(exact,3));

%% HELMHOLTZ WITH K = 2  AND NEUMANN BC

K = 2;

% Example 20:
exact = ballfun(@(r, lam, th)1,'polar');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 0;  % Neumann condition
u = helmholtz_neumann(f, K, bc1, 38, 17, 24);
pass(20) = isequal(u,exact);

% Example 21:
exact = ballfun(@(r, lam, th)r.^2.*sin(th).^2,'polar');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 2*sin(th).^2;  % Neumann condition
u = helmholtz_neumann(f, K, bc1,38,17,24);
pass(21) = isequal(u,exact);

% Example 22:
exact = ballfun(@(r, lam, th)r.^4.*sin(th).^2,'polar');
f = laplacian(exact) + K^2*exact;
bc1 = @(lam, th) 4*sin(th).^2;  % Neumann condition
u = helmholtz_neumann(f, K, bc1,38,17,24);
pass(22) = isequal(u,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end
