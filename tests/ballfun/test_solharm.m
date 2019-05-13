function pass = test_solharm( pref ) 
% Test the equality for P^m_l, 0 <= l <= n and -m <= l <= m

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e6*pref.techPrefs.chebfuneps;

n = 10;
Max_difference = 0;
Max_norm = 0;
Max_laplacian = 0;

for l = 0:n
    for m = -l:l
        Yml = ballfun.solharm(l,m);
        Y = Yml(1,:,:,'spherical')/sqrt(2*l+3);
        Z = spherefun.sphharm(l,m);
        % Check that the Legendre polynomials is the same as in Spherefun
        Max_difference = max(norm(Y-Z),Max_difference);
        % Check that the 2-norm is 1
        Max_norm = max(abs(norm(Yml)-1),Max_norm);
        % Check that they are solutions to the Laplace equation
        Max_laplacian = max(norm(laplacian(Yml)),Max_laplacian);
    end
end
pass(1) = Max_difference < tol;
pass(2) = Max_norm < tol;
pass(3) = Max_laplacian < tol;

% Example 1 : Y^0_0
m = 0; l = 0;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)0.5*sqrt(1/pi)*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(4) = norm(Exact - f) < tol;

% Example 2 : Y^-1_1
m = -1; l = 1;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)0.5*sqrt(3/(2*pi))*r.^l.*exp(1i*m*lam).*sin(th)*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(5) = norm(Exact - f) < tol;

% Example 4 : Y^0_1
m = 0; l = 1;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)0.5*sqrt(3/pi)*r.^l.*exp(1i*m*lam).*cos(th)*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(6) = norm(Exact - f) < tol;

% Example 5 : Y^1_1
m = 1; l = 1;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)-0.5*sqrt(3/(2*pi))*r.^l.*exp(1i*m*lam).*sin(th)*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(7) = norm(Exact - f) < tol;

% Example 6 : Y^-2_2
m = -2; l = 2;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)0.25*sqrt(15/(2*pi))*r.^l.*exp(1i*m*lam).*sin(th).^2*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(8) = norm(Exact - f) < tol;

% Example 7 : Y^-1_2
m = -1; l = 2;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)0.5*sqrt(15/(2*pi))*r.^l.*exp(1i*m*lam).*sin(th).*cos(th)*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(9) = norm(Exact - f) < tol;

% Example 8 : Y^0_2
m = 0; l = 2;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)0.25*sqrt(5/pi)*r.^l.*exp(1i*m*lam).*(3*cos(th).^2-1)*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(10) = norm(Exact - f) < tol;

% Example 9 : Y^1_2
m = 1; l = 2;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)-0.5*sqrt(15/(2*pi))*r.^l.*exp(1i*m*lam).*sin(th).*cos(th)*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(11) = norm(Exact - f) < tol;

% Example 10 : Y^2_2
m = 2; l = 2;
Normalization = sqrt(2*l+3);
Exact = ballfun(@(r,lam,th)0.25*sqrt(15/(2*pi))*r.^l.*exp(1i*m*lam).*sin(th).^2*Normalization, 'spherical');
f = ballfun.solharm(l,m,'complex');
pass(12) = norm(Exact - f) < tol;

% Example 11: real solid harmonics
f = ballfun.solharm(3,2);
pass(13) = abs(norm(f)-1) < tol;

% Example 12: Y^0_5
f = ballfun.solharm(5,0);
pass(14) = abs(norm(f)-1) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end