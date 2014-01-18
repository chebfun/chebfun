function [bol, wzero] = singleSignTest(F)
%SINGLESIGNTEST(F) heuristic check to see if F changes sign.
% 
% SINGLESIGNTEST(F) returns 1 if the values of F on a tensor grid are of the
% same sign. 
%
% [BOL, WZERO] = SINGLESIGNTEST(F), WZERO = 1 if a zero has been found. 
%
% See also ABS, SQRT, LOG. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

bol = 0 ;   % assume false. 
wzero = 0;   % assume no zeros.

X = chebpolyval2(F);  % evaluate on a grid use FFTs. 

if ( all(all( X >=0 )))  % If all values are nonnegative 
    bol = 1;  
elseif ( all(all( X <=0 )) )% If all values are not positive
    bol = 1; 
end

if ( any ( any ( X == 0 ) ) ) 
    wzero = 1; 
end

end