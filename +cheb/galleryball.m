function f = galleryball(name,varargin)
%CHEB.GALLERYBALL   Ballfun example functions.
%   CHEB.GALLERYBALL(NAME) returns a ballfun function corresponding to
%   NAME.  See the listing below for available names.
%
%   gaussian   Gaussian function on the ball centered at (-0.5,0,0)
%   solharm    Solid harmonics of degree 5 and order 3

% If the user did not supply an input, return a function chosen at random
% from the gallery.
if ( nargin == 0 )
    names = {'gaussian', 'solharm'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch lower(name)
    
    % Gaussian function centered at (-0.5,0,0)
    case 'gaussian'
        fa = @(x,y,z,xc,yc,zc) exp(-20*((x+0.5).^2 + y.^2 + z.^2));
        f = ballfun(@(x,y,z) fa(x,y,z));

    % Solid harmonics function
    case 'solharm'
        f = ballfun.solharm(5,3);
         
     % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALL:unknown:unknownFunction', ...
            'Unknown function.')
end
end
