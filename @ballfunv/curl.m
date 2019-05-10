function G = curl(F)
%CURL Curl of a BALLFUNV.
%   G = CURL(F) returns the BALLFUNV of the curl of F.
%
% See also DIV. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( F )
    G = ballfunv();
    return
end

% Extract the components of a BALLFUNV:
Fc = F.comp;

% Formula for 3D curl
G = [ diff(Fc{3},2) - diff(Fc{2},3); ...
      diff(Fc{1},3) - diff(Fc{3},1); ...
      diff(Fc{2},1) - diff(Fc{1},2) ];
end