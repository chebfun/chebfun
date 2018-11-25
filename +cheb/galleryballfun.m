function f = galleryballfun(name,varargin)
%CHEB.GALLERYBALLFUN   Ballfun example functions.
%   CHEB.GALLERYBALLFUN(NAME) returns a ballfun function corresponding to
%   NAME.  See the listing below for available names.
%
%   zero          Zero function

% The main switch statement.
switch lower(name)

    % Zero function
    case 'zero'
        f = ballfun(0,'coeffs');
         
     % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALLFUN:unknown:unknownFunction', ...
            'Unknown function.')
end
end
