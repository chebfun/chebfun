function varargout = chebfun3eps(val)
%CHEBFUN3EPS   Set the default value of the chebfun3eps preference.
%   CHEBFUN3EPS VAL, or CHEBFUN3EPS(VAL) sets the default value of the
%   chebfun3eps preference to the specified value VAL. CHEBFUN3EPS(VAL) is 
%   equivalent to CHEBFUNPREF.SETDEFAULTS({'cheb3Prefs', 'chebfun3eps'}, VAL).
%
%   CHEBFUN3EPS factory, or CHEBFUN3EPS('factory') sets the default 
%   chebfun3eps preference to the factory value. CHEBFUN3EPS('factory') is
%   equivalent to CHEBFUNPREF.SETDEFAULTS({'cheb3prefs', 'chebfun3eps'},
%   'factory').
%
%   CHEBFUN3EPS prints the current default value of the chebfun3eps 
%   preference.
%
%   If changing the preference is needed only for a single construction, 
%   calling constructor with the 'eps' flag is a better option.
%
% See also CHEBFUNEPS and CHEBFUN2EPS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )
    % Return current chebfun3eps:
    varargout{1} = chebfunpref().cheb3Prefs.chebfun3eps;
else
    if ( strcmpi(val, 'factory') )
        chebfunpref.setDefaults({'cheb3Prefs','chebfun3eps'}, 'factory');
    else
        eval(['val=',val,';'])
        if ( isnumeric(val) )
            chebfunpref.setDefaults({'cheb3Prefs','chebfun3eps'}, val);
        else
            error('CHEBFUN:chebfun3eps:unknownOption',...
                'Unknown CHEBFUN3EPS option.')
        end
    end
end

end
