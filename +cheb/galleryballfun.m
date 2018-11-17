function f = galleryballfun(name,varargin)
%CHEB.GALLERYBALLFUN   Ballfun example functions.
%   CHEB.GALLERYBALLFUN(NAME) returns a ballfun function corresponding to
%   NAME.  See the listing below for available names.
%
%   random        Interpolant through random data in CFF grid
%   zero          Zero function

% The main switch statement.
switch lower(name)
    
    % Interpolant through random data
    case 'random'
        f = randnfunball(5);
        if nargin > 1
            S = varargin{1};
            F = coeffs3(f,S(1),S(2),S(3));
            f = ballfun(F,'coeffs');
        end

    % Zero function
    case 'zero'
        f = ballfun(@(r,lam,th)0);
         
     % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERYBALLFUN:unknown:unknownFunction', ...
            'Unknown function.')
end
end
