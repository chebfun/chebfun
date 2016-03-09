function varargout = chebfun2eps(val)
%CHEBFUN2EPS   Set chebfun2eps globally.
%   CHEBFUN2EPS val, or CHEBFUN2EPS(val) sets the value of chebfun2eps to 
%   be equal to the specified value val.
%
%   CHEBFUN2EPS factory, or CHEBFUN2EPS('factory') sets the chebfun2eps to 
%   be equal to the factory value.
%
%   CHEBFUN2EPS prints the current value of chebfun2eps.
%
% See also CHEBFUNEPS and CHEBFUN3EPS.

if ( nargin == 0 )
    % Return current chebfun2eps:
    chebfun2epsVal = chebfunpref().cheb2Prefs.chebfun2eps;
    varargout{1} = chebfun2epsVal;
    
else
    
    if ( strcmpi(val, 'factory') )
        prefs = chebfunpref.getFactoryDefaults();
        val = prefs.cheb2Prefs.chebfun2eps;
        chebfunpref.setDefaults({'cheb2Prefs','chebfun2eps'},val);
        
    elseif double(val)
        if ischar(val)
            val = str2double(val);
        end
        chebfunpref.setDefaults({'cheb2Prefs','chebfun2eps'},val);
        
    else
        error('CHEBFUN:CHEBFUB2:chebfun2eps:UnknownOption',...
            'Unknown chebfun2eps option.')
    end
end

end