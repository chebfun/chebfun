function f = randnfundisk(dt)
%RANDNFUNDISK   Random smooth function on the unit disk
%   F = RANDNFUNDISK(DT) returns a smooth DISKFUN of maximum
%   frequency about 2pi/DT and standard normal distribution N(0,1)
%   at each point.  F is obtained by restricting the output of
%   RANDNFUN2 to the unit disk.
%
%   RANDNFUNDISK() uses the default value DT = 1.
%
% Examples:
%
%   f = randnfundisk(0.5); plot(f)
%   colormap([0 0 0; 1 0 0]), caxis(norm(caxis,inf)*[-1 1]), axis off
%
% See also RANDNFUN, RANDNFUN2, RANDFUNSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 0
    dt = 1;
end
 
fsquare = randnfun2(dt);
f = diskfun(@(x,y) fsquare(x,y));

end
