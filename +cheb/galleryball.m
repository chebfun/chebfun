function f = galleryball(name,varargin)
%CHEB.GALLERYBALL   Ballfun example functions.
%   CHEB.GALLERYBALL(NAME) returns a ballfun function corresponding to
%   NAME.  See the listing below for available names.
%
%   gaussian    Gaussian function on the ball centered at (-0.5,0,0)
%   peaks       Peaks like function on the ball taken from the geopeaks
%               function in the MATLAB mapping toolbox.
%   solharm     Solid harmonics of degree 5 and order 3

% If the user did not supply an input, return a function chosen at random
% from the gallery.
if ( nargin == 0 )
    names = {'gaussian', 'peaks', 'solharm'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch lower(name)
    
    % Gaussian function centered at (-0.5,0,0)
    case 'gaussian'
        fa = @(x,y,z,xc,yc,zc) exp(-20*((x+0.5).^2 + y.^2 + z.^2));
        f = ballfun(@(x,y,z) fa(x,y,z));
    
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
         
     % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALL:unknown:unknownFunction', ...
            'Unknown function.')
end
end
