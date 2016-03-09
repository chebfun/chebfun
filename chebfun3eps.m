function varargout = chebfun3eps(val)
%CHEBFUN3EPS   Set chebfun3eps globally.
%   CHEBFUN3EPS val, or CHEBFUN2EPS(val) sets the value of chebfun3eps to 
%   the specified value val.
%
%   CHEBFUN3EPS('factory') sets the chebfun3eps to be equal to 
%   the factory value val.
%
%   CHEBFUN3EPS prints the current value of chebfun3eps.
%
% See also CHEBFUNEPS and CHEBFUN2EPS.

if ( nargin == 0 )
    % Return current chebfun3eps:
    chebfun3epsVal = chebfunpref().cheb3Prefs.chebfun3eps;
    varargout{1} = chebfun3epsVal;
    
else
    
    if ( strcmpi(val, 'factory') )
        prefs = chebfunpref.getFactoryDefaults();
        val = prefs.cheb3Prefs.chebfun3eps;
        chebfunpref.setDefaults({'cheb3Prefs','chebfun3eps'},val);
        
    elseif double(val)
        if ischar(val)
            val = str2double(val);
        end
        chebfunpref.setDefaults({'cheb3Prefs','chebfun3eps'},val);
        
    else
        error('CHEBFUN:CHEBFUB3:chebfun3eps:UnknownOption',...
            'Unknown chebfun3eps option.')
    end
end

end