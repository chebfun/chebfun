function varargout = fevalm(f, r, lambda, theta)
%FEVAL   Evaluate a BALLFUNV
%   FEVAL(F, R, L, T) evaluates a BALLFUNV F at the points (R,L,T) in spherical coordinates
%   at a tensor-product grid R x L x T.
%
% See also SUBSREF. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f )
    varargout = {};
    return
end

F = f.comp;

% Evaluate at the tensor grid
valsX = fevalm(F{1}, r, lambda, theta);
valsY = fevalm(F{2}, r, lambda, theta);
valsZ = fevalm(F{3}, r, lambda, theta);

varargout = {valsX, valsY, valsZ};
end
