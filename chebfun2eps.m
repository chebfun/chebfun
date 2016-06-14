function varargout = chebfun2eps(val)
%CHEBFUN2EPS   Set the default value of the chebfun2eps preference.
%   CHEBFUN2EPS VAL, or CHEBFUN2EPS(VAL) sets the default value of the
%   chebfun2eps preference to the specified value VAL. CHEBFUN2EPS(VAL) is 
%   equivalent to CHEBFUNPREF.SETDEFAULTS({'cheb2Prefs', 'chebfun2eps'}, VAL).
%
%   CHEBFUN2EPS factory, or CHEBFUN2EPS('factory') sets the default 
%   chebfun2eps preference to the factory value. CHEBFUN2EPS('factory') is
%   equivalent to CHEBFUNPREF.SETDEFAULTS({'cheb2prefs', 'chebfun2eps'},
%   'factory').
%
%   CHEBFUN2EPS prints the current default value of the chebfun2eps 
%   preference.
%
%   If changing the preference is needed only for a single construction, 
%   calling constructor with the 'eps' flag is a better option.
%
% See also CHEBFUNEPS and CHEBFUN3EPS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )
    % Return current chebfun2eps:
    varargout{1} = chebfunpref().cheb2Prefs.chebfun2eps;
else
    if ( strcmpi(val, 'factory') )
        chebfunpref.setDefaults({'cheb2Prefs','chebfun2eps'}, 'factory');
    else
        eval(['val=',val,';'])
        if ( isnumeric(val) )
            chebfunpref.setDefaults({'cheb2Prefs','chebfun2eps'}, val);
        else
            error('CHEBFUN:chebfun2eps:unknownOption',...
                'Unknown CHEBFUN2EPS option.')
        end
    end
end


end
