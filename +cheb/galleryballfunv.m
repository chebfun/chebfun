function v = galleryballfunv(name,S)
%CHEB.GALLERYBALLFUNV   Ballfunv example.
%   CHEB.GALLERYBALLFUNV(NAME, S) returns a ballfunv corresponding to
%   NAME with size S.  See the listing below for available names.
%
%   random         Smooth random vector field
%   zero           Zero vector field

% The main switch statement.
switch lower(name)
    
    case 'random'
        v = randnfunballv(5,S);
        
    % Zero vector field
    case 'zero'
        f = cheb.galleryballfun('zero',S);
        v = ballfunv(f,f,f);
        
     % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALLFUNV:unknown:unknownFunction', ...
            'Unknown function.')
end
end


% Compute the initial vorticity for the NS equations
function W = ns_vorticity(S)
% Two dipole cylindrical spinning tops that will collide against the
% spherical walls in a NS simulation:
x1 = .4; 
x2 = -.4;
% Change this, not correct as written:
Wx = @(x, y, z) 2*(exp(-30*((x-x1).^2+y.^2)) - exp(-30*((x-x2).^2+y.^2))).*exp(-30*z.^2);
Wy = 0; 
Wz = 0;

% Convert to spherical coordinates: 
x = @(r, th, lam) r.*sin(th).*sin(lam);
y = @(r, th, lam) r.*sin(th).*cos(lam);
z = @(r, th, lam) r.*cos(th);

% Vector field in spherical coordinates: 
Wr = ballfun( @(r, th, lam) sin(th).*cos(lam).*Wx(x(r,th,lam),y(r,th,lam),z(r,th,lam)), S );
Wth = ballfun( @(r, th, lam) cos(th).*cos(lam).*Wx(x(r,th,lam),y(r,th,lam),z(r,th,lam)), S );
Wlam = ballfun( @(r, th, lam) -sin(lam).*Wx(x(r,th,lam),y(r,th,lam),z(r,th,lam)), S );

% Make a ballfunv object:
P = ballfunv( Wr , Wth , Wlam );
W = -laplace(P);
end

% 2nd initial vorticity
function W = ns_vorticity_2(S)
    alpha = 5.763459196894548;
    lambda_0 = 10.53345116848308;
    % psi = @(r,lam,th)r^2*sin(th)^2*(lambda_0/(alpha^2)+1/(r^(3/2))*besselj(3/2,alpha*r));
    
    % Derivative of the bessel function
    % d/dz besselj(nu,z) = (nu*besselj(nu, z))/z - besselj(nu + 1, z)
    
    % cylindrical to spherical
    % rho, lambda, z -> r, lambda, theta
    % r = sqrt(rho^2+z^2)
    % theta = tan^-1(rho/z)
    % lambda = lambda
    
    % dpsi/dR
    % Psi_R = @(r,lam,th)sin(th)^2*(2*r*lambda_0/(alpha^2)+besselj(3/2,alpha*r)/(2*sqrt(r))+...
    %                   alpha*sqrt(r)*((3/2*besselj(3/2, alpha*r))/(alpha*r) - besselj(3/2 + 1, alpha*r)));
    % dpsi/dth
    % Psi_th = @(r,lam,th)2*r^2*sin(th)*cos(th)*(lambda_0/(alpha^2)+1/(r^(3/2))*besselj(3/2,alpha*r));
    
    % 1/r*dR/dr*dpsi/dR
    % a = @(r,lam,th)sin(th)^2*(2*lambda_0/(alpha^2)+besselj(3/2,alpha*r)/(2*r^(3/2))+...
    %                    3/2*r^(-3/2)*besselj(3/2, alpha*r) - alpha*besselj(3/2 + 1, alpha*r)/sqrt(r));
    % 1/r*dR/dz*dpsi/dR
    % b = @(r,lam,th)cos(th)*sin(th)*(2*lambda_0/(alpha^2)+besselj(3/2,alpha*r)/(2*r^(3/2))+...
    %                   3/2*r^(-3/2)*besselj(3/2, alpha*r) - alpha*besselj(3/2 + 1, alpha*r)/sqrt(r));
   
    % 1/r*dth/dr*dpsi/dth
    % c = @(r,lam,th)2*cos(th)^2*(lambda_0/(alpha^2)+1/(r^(3/2))*besselj(3/2,alpha*r));
    
    % 1/r*dth/dz*dpsi/dth
    % d = @(r,lam,th)-2*sin(th)*cos(th)*(lambda_0/(alpha^2)+1/(r^(3/2))*besselj(3/2,alpha*r));
    
    % Vector field in cylindrical coordinates
    % -(b+d)
    Vcyl1 = @(r,lam,th)sin(th)*cos(th)*alpha*besselj(3/2+1,alpha*r)/sqrt(r);
    Vcyl2 = @(r,lam,th)alpha*sin(th)*(r*lambda_0/(alpha^2)+1/sqrt(r)*besselj(3/2,alpha*r));
    % a+c
    Vcyl3 = @(r,lam,th)2*lambda_0/(alpha^2)+2*r^(-3/2)*besselj(3/2,alpha*r)-alpha*sin(th)^2*besselj(3/2+1,alpha*r)/sqrt(r);
    
    % Transform a vector field in cylindrical coordinates to spherical
    % coordinates

%     Mcartcyl =
%  
% [  cos(th), sin(th), 0]
% [ -sin(th), cos(th), 0]
% [        0,       0, 1]
% Mcartsph =
%  
% [ cos(lam)*sin(th), sin(lam)*sin(th),  cos(th)]
% [ cos(lam)*cos(th), cos(th)*sin(lam), -sin(th)]
% [        -sin(lam),         cos(lam),        0]

% M = Mcartsph*Mcartcyl^-1 = Mcartsph*Mcartcyl.'
%  
% M =
%  
% [ sin(lam)*sin(th)^2 + cos(lam)*cos(th)*sin(th), cos(th)*sin(lam)*sin(th) - cos(lam)*sin(th)^2,  cos(th)]
% [ cos(lam)*cos(th)^2 + sin(lam)*sin(th)*cos(th), cos(th)^2*sin(lam) - cos(lam)*cos(th)*sin(th), -sin(th)]
% [           cos(lam)*sin(th) - cos(th)*sin(lam),           cos(lam)*cos(th) + sin(lam)*sin(th),        0]
 
    Vsph1 = @(r,lam,th)(sin(lam)*sin(th)^2+cos(lam)*cos(th)*sin(th))*Vcyl1(r,lam,th)...
                     + (sin(lam)*cos(th)*sin(th)-cos(lam)*sin(th)^2)*Vcyl2(r,lam,th)...
                     + cos(th)*Vcyl3(r,lam,th);
    Vsph2 = @(r,lam,th)(cos(lam)*cos(th)^2+sin(lam)*cos(th)*sin(th))*Vcyl1(r,lam,th)...
                     + (sin(lam)*cos(th)^2-cos(lam)*cos(th)*sin(th))*Vcyl2(r,lam,th)...
                     - sin(th)*Vcyl3(r,lam,th);
    Vsph3 = @(r,lam,th)(cos(lam)*sin(th)-cos(th)*sin(lam))*Vcyl1(r,lam,th)...
                     + (cos(lam)*cos(th)+sin(th)*sin(lam))*Vcyl2(r,lam,th);
    % Create ballfun and ballfunv
    Vsph1 = ballfun(Vsph1,S);
    Vsph2 = ballfun(Vsph2,S);
    Vsph3 = ballfun(Vsph3,S);
    V = ballfunv(Vsph1,Vsph2,Vsph3);
    W = curl(V);
end
