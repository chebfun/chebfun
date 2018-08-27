function f = galleryballfun(name,varargin)
%CHEB.GALLERYBALLFUN   Ballfun example functions.
%   CHEB.GALLERYBALLFUN(NAME, S) returns a ballfun function corresponding to
%   NAME with size S.  See the listing below for available names.
%
%   random        Interpolant through random data in CFF grid
%   zero          Zero function

% The main switch statement.
switch lower(name)
    
    % Interpolant through random data
    case 'random'
        f = randnfunball(5,varargin{1});

    % Zero function
    case 'zero'
        f = ballfun(@(r,lam,th)0);
         
     % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALLFUN:unknown:unknownFunction', ...
            'Unknown function.')
end
end
