function f = galleryball(name,varargin)
%CHEB.GALLERYBALL   Ballfun example functions.
%   CHEB.GALLERYBALL(NAME) returns a ballfun function corresponding to
%   NAME.  See the listing below for available names.
%
%   deathstar   A function resembling the Death Star.
%   gaussian    Gaussian function on the ball centered at (-0.5,0,0)
%   moire       Moire pattern from waves generated at two point sources.
%   peaks       Peaks like function on the ball taken from the geopeaks
%               function in the MATLAB mapping toolbox.
%   solharm     Solid harmonics of degree 5 and order 3
%   stripes     Alternating striped pattern.



% If the user did not supply an input, return a function chosen at random
% from the gallery.
if ( nargin == 0 )
    names = {'deathstar','gaussian','moire','peaks','solharm','stripes'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch lower(name)
    
    % A function resembling the Death Star
    case 'deathstar'
        fa = @(x,y,z) -(exp(-30*((y+sqrt(3)/2).^2 + x.^2 + (z-1/2).^2)) + exp(-100*z.^2));
        f = ballfun(fa);
        
    % Gaussian function centered at (-0.5,0,0)
    case 'gaussian'
        fa = @(x,y,z,xc,yc,zc) exp(-20*((x+0.5).^2 + y.^2 + z.^2));
        f = ballfun(@(x,y,z) fa(x,y,z));
    
    % An interesting Moire pattern generated by two sources on the ball
    case 'moire'
        % Centers of the beacons
        boise = [-116.237651 43.613739]*pi/180;
        oxford = [-1.257778 51.751944]*pi/180;
        % ithaca = [-76.5 42.443333]*pi/180;
        % stellenbosh = [18.86 -33.92]*pi/180;
        [xb,yb,zb] = sph2cart(boise(1),boise(2),1);
        [xo,yo,zo] = sph2cart(oxford(1),oxford(2),1);
        % Pick the number of oscillations and make each of the "waves" 
        % vanish at the anti-podal points from their centers.
        omega = besselroots(0,30); omega = omega(end)/2;
        % Use a combination of the J0 bessel functions centered at Boise
        % and Oxford to generate the Moire pattern.
        fa = @(x,y,z,omega) 2 + besselj(0,omega*sqrt((x-xb).^2+(y-yb).^2+(z-zb).^2)) + ...
            2 + besselj(0,omega*sqrt((x-xo).^2+(y-yo).^2+(z-zo).^2));
        f = ballfun(@(x,y,z) fa(x,y,z,omega));
        
    % A "peaks-like" function for the ball
    case 'peaks'
        fa = @(x,y,z) 8*(1-x).^2.*exp(-4*(x - 0.059).^2 - 2*(y + 0.337).^2 - 2*(z + 0.940).^2) - ...
            30*(z/10 - x.^3 - y.^5) .* exp(-3*(x - 0.250).^2 - 2*(y - 0.433).^2 - 3*(z - 0.866).^2) + ...
            (20*y - 8*z.^3) .* exp(-2*(x + 0.696).^2 - 3*(y + 0.123).^2 - 2*(z - 0.707).^2) + ...
            (7*y - 10*x + 10*z.^3) .* exp(-3*(x - 0.296).^2 - 3*(y + 0.814).^2 - 3*(z + 0.5).^2);
        f = ballfun(fa);
        
    % Solid harmonics function
    case 'solharm'
        f = ballfun.solharm(5,3);
        
    % A zipper-like stripe pattern for the sphere    
    case 'stripes'
        fa = @(x,y,z) (1 + cos(10*pi*x)).*exp(-exp(-20*z)) + (1 - cos(10*pi*x)).*exp(-exp(20*z));
        f = ballfun(fa);
        
        
    % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALL:unknown:unknownFunction', ...
            'Unknown function.')
end
end