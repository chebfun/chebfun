function v = mean3(f)
%MEAN3   Mean of a CHEBFUN3.
%   V = MEAN3(F) returns the mean of a CHEBFUN3 object F, i.e., sum3(F) / V,
% 	where V is the volume of the domain of F.
%
% See also CHEBFUN3/MEAN, CHEBFUN3/MEAN2 and CHEBFUN3/STD3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) ) 
    v = [];
    return
end

% Apply the formula: 
v = sum3(f) / domainvolume(f);  

end