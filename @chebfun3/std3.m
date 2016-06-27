function v = std3(f)
%STD3   Standard deviation of a CHEBFUN3.
%   V = STD3(F) computes the standard deviation of the CHEBFUN3 object F, 
%   i.e.,
%             STD3(F)^2 = 1 / V*sum3(|F(x,y,z) - M|^2)
%
%   where V is the volume of the domain of F and M is the mean of F.
%
% See also CHEBFUN3/MEAN, CHEBFUN3/MEAN2, CHEBFUN3/MEAN3, CHEBFUN3/STD and 
% CHEBFUN3/STD2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = f - mean3(f);
v = sqrt(mean3(h .* conj(h)));

end