function varargout = chebfuneps(val)
%CHEBFUNEPS   Set the default value of the chebfuneps preference.
%   CHEBFUNEPS VAL, or CHEBFUNEPS(VAL) sets the value of chebfuneps to be 
%   equal to the specified value VAL. CHEBFUNEPS(VAL) is equivalent to 
%   CHEBFUNPREF.SETDEFAULTS('chebfuneps', VAL).
%
%   CHEBFUNEPS factory, or CHEBFUNEPS('factory') sets the chebfuneps to be equal to 
%   the factory value.
%
%   CHEBFUNEPS prints the current value of chebfuneps.
%
%   If changing the preference is needed only for a single construction, 
%   calling constructor with the 'eps' flag is a better option.
%
% See also CHEBFUN2EPS.

if ( nargin == 0 )
    % Return current chebfuneps:
    chebfunepsVal = chebfunpref().techPrefs.chebfuneps
    varargout{1} = chebfunepsVal;
    
else
    
    if ( strcmpi(val, 'factory') )
        prefs = chebfunpref.getFactoryDefaults();
        val = prefs.techPrefs.chebfuneps;
        chebfunpref.setDefaults('chebfuneps',val);
        
    elseif isnumeric(val)
        if ischar(val)
            val = str2double(val);
        end
        chebfunpref.setDefaults('chebfuneps',val);
        
    else
        error('CHEBFUN:chebfuneps:UnknownOption',...
            'Unknown CHEBFUNEPS option.')
    end
end

end