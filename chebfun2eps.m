function varargout = chebfun2eps(val)
%CHEBFUN2EPS   Set the default value of the chebfun2eps preference.
%   CHEBFUN2EPS VAL, or CHEBFUN2EPS(VAL) sets the value of chebfun2eps to 
%   be equal to the specified value VAL. CHEBFUN2EPS(VAL) is equivalent to 
%   CHEBFUNPREF.SETDEFAULTS({'cheb2prefs','chebfun2eps'}, VAL).
%
%   CHEBFUN2EPS factory, or CHEBFUN2EPS('factory') sets the chebfun2eps to 
%   be equal to the factory value.
%
%   CHEBFUN2EPS prints the current value of chebfun2eps.
%
%   If changing the preference is needed only for a single construction, 
%   calling constructor with the 'eps' flag is a better option.
%
% See also CHEBFUNEPS.

if ( nargin == 0 )
    % Return current chebfun2eps:
    chebfun2epsVal = chebfunpref().cheb2Prefs.chebfun2eps;
    varargout{1} = chebfun2epsVal;
    
else

    if ( strcmpi(val, 'factory') )
        prefs = chebfunpref.getFactoryDefaults();
        val = prefs.cheb2Prefs.chebfun2eps;
        chebfunpref.setDefaults({'cheb2Prefs','chebfun2eps'},val);
        
    elseif isnumeric(val)
        if ischar(val)
            val = str2double(val);
        end
        chebfunpref.setDefaults({'cheb2Prefs','chebfun2eps'},val);
        
    else
        error('CHEBFUN:chebfun2eps:UnknownOption',...
            'Unknown CHEBFUN2EPS option.')
    end
end

end