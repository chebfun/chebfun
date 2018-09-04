function F = curl(G)
%CURL Curl of a BALLFUNV in cartesian coordinates.
%   CURL(G) is the curl of the BALLFUNV w.
%
% See also DIV. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Gc = G.comp;

F = [ diff(Gc{3},2,"cart") - diff(Gc{2},3,"cart"); ...
      diff(Gc{1},3,"cart") - diff(Gc{3},1,"cart"); ...
      diff(Gc{2},1,"cart") - diff(Gc{1},2,"cart") ];
end