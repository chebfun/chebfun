function iscart = isCartesian( varargin )
% COORDSETTING  This function searches for the 'polar' flag, which 
% may be present in a call to the diskfun constructor or in feval. 
% When the 'polar' flag is present, evaluations should be done
% with respect to polar coordinates instead of Cartesian coordinates. 
%
% ISCART = COORDSETTING(varargin). If the flag 'polar' is included
% in varargin, ISCART = 0. Otherwise, ISCART = 1. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


iscart = 1; 
% search for user-supplied 'polar' flag in arguments: 
isPolar = find(strcmp(varargin, 'polar'));
if ( any( isPolar ) )
    iscart = 0; 
end


end




