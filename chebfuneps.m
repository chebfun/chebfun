function varargout = chebfuneps(val)
%CHEBFUNEPS   Set chebfuneps globally.
%   CHEBFUNEPS val, or CHEBFUNEPS(val) sets the value of chebfuneps to be equal to the 
%   specified value val.
%
%   CHEBFUNEPS factory, or CHEBFUNEPS('factory') sets the chebfuneps to be equal to 
%   the factory value.
%
%   CHEBFUNEPS prints the current value of chebfuneps.
%
% See also CHEBFUN2EPS and CHEBFUN3EPS.

if ( nargin == 0 )
    % Return current chebfuneps:
    chebfunepsVal = chebfunpref().techPrefs.chebfuneps
    varargout{1} = chebfunepsVal;
    
else
    
    if ( strcmpi(val, 'factory') )
        prefs = chebfunpref.getFactoryDefaults();
        val = prefs.techPrefs.chebfuneps;
        chebfunpref.setDefaults('chebfuneps',val);
        
    elseif double(val)
        if ischar(val)
            val = str2double(val);
        end
        chebfunpref.setDefaults('chebfuneps',val);
        
    else
        error('CHEBFUN:chebfuneps:UnknownOption',...
            'Unknown chebfuneps option.')
    end
end

end