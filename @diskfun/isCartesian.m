function iscart = isCartesian( varargin )
% COORDSETTING     Decide if the function is in polar or Cartesian:
%   ISCART = COORDSETTING(F ) returns ISCART = 1 unless the flag 'polar'
% is present. Then, ISCART=0.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


iscart = 1; 
% search for user-supplied 'polar' flag in arguments: 
isPolar = find(strcmp(varargin, 'polar'));
if ( any( isPolar ) )
    iscart = 0; 
end


end




