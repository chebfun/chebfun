function vals = feval(varargin)
%FEVAL   Evaluate a BALLFUNV
%   FEVAL(F, R, L, T) evaluates a BALLFUNV F at the points (R,L,T) in spherical coordinates
%   at a tensor-product grid R x L x T.
%
% See also SUBSREF. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = varargin{1};

F = f.comp;

V1 = feval(F{1}, varargin{2:end});
V2 = feval(F{2}, varargin{2:end});
V3 = feval(F{3}, varargin{2:end});
vals = zeros(size(V1,3));
vals(:,:,:,1) = V1;
vals(:,:,:,2) = V2;
vals(:,:,:,3) = V3;
end
