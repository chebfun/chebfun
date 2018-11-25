function v = galleryballfunv(name)
%CHEB.GALLERYBALLFUNV   Ballfunv example.
%   CHEB.GALLERYBALLFUNV(NAME) returns a ballfunv corresponding to
%   NAME with size S.  See the listing below for available names.
%
%   zero            Zero vector field

% The main switch statement.
switch lower(name)
    
    % Zero vector field
    case 'zero'
        f = cheb.galleryballfun('zero');
        v = ballfunv(f,f,f);
        
     % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALLFUNV:unknown:unknownFunction', ...
            'Unknown function.')
end
end
