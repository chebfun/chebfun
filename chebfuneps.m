function varargout = chebfuneps(val)
%CHEBFUNEPS   Set the default value of the chebfuneps preference.
%   CHEBFUNEPS VAL, or CHEBFUNEPS(VAL) sets the default value of the 
%   chebfuneps preference to the specified value VAL. CHEBFUNEPS(VAL) is 
%   equivalent to CHEBFUNPREF.SETDEFAULTS('chebfuneps', VAL).
%
%   CHEBFUNEPS factory, or CHEBFUNEPS('factory') sets the default 
%   chebfuneps preference to the factory value. This is equivalent to 
%   chebfunpref.setDefaults('chebfuneps', 'factory').
%
%   CHEBFUNEPS prints the current default value of the chebfuneps 
%   preference.
%
%   If changing the preference is needed only for a single construction, 
%   calling constructor with the 'eps' flag is a better option.
%
% See also CHEBFUN2EPS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )
    % Return current chebfuneps:
    varargout{1} = chebfunpref().techPrefs.chebfuneps;
else
    if ( strcmpi(val, 'factory') )
        chebfunpref.setDefaults('chebfuneps', 'factory');
    else
        eval(['val=',val,';'])
        if ( isnumeric(val) )
            chebfunpref.setDefaults('chebfuneps', val);
        else
            error('CHEBFUN:chebfuneps:unknownOption',...
                'Unknown CHEBFUNEPS option.')
        end
    end
end

end
