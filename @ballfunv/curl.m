function G = curl(F)
%CURL Curl of a BALLFUNV.
%   S = CURL(F) returns the DISKFUN of the curl of F.
%
% See also DIV. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract the components of a BALLFUNV:
Fc = F.comp;

% Formula for 3D curl
G = [ diff(Fc{3},2) - diff(Fc{2},3); ...
      diff(Fc{1},3) - diff(Fc{3},1); ...
      diff(Fc{2},1) - diff(Fc{1},2) ];
end