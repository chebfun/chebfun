function f = curl(F)
%CURL  curl of a chebfun2v
% 
% S = CURL(F) returns the chebfun2 of the curl of F. If F is a chebfun2v 
% with two components then it returns the chebfun2 representing 
%
%         CURL(F) = F(2)_x - F(1)_y,
%
% where F = (F(1),F(2)).  If F is a chebfun2v with three components then it
% returns the chebfun2v representing the 3D curl operation. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.  

components = F.components; 
if ( F.nComponents == 2 ) % do the curl of a 2-vector. 
    f = diff(components(2),1,2) - diff(components(1),1,1);
else   % do the curl of a 3-vector.
    curlVector = [diff(component(3),1,1)...
        -diff(components(3),1,2)...
        diff(components(2),1,2) - diff(components(1),1,1)];
    
    f.components = curlVector;  
end

end