function f = galleryballfun(name,varargin)
%CHEB.GALLERYBALLFUN   Ballfun example functions.
%   CHEB.GALLERYBALLFUN(NAME) returns a ballfun function corresponding to
%   NAME.  See the listing below for available names.
%
%   solharm    Solid harmonics of degree 5 and order 3

% If the user did not supply an input, return a function chosen at random
% from the gallery.
if ( nargin == 0 )
    names = {'solharm'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch lower(name)

    % Solid harmonics function
    case 'solharm'
        f = ballfun.solharm(5,3);
         
     % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALLFUN:unknown:unknownFunction', ...
            'Unknown function.')
end
end
