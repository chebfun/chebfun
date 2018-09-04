function vals = feval(f, r, lambda, theta)
%FEVAL   Evaluate a BALLFUNV
%   FEVAL(F, R, L, T) evaluates a BALLFUNV F at the points (R,L,T) in spherical coordinates
%   at a tensor-product grid R x L x T.
%
% See also SUBSREF. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;

% Get the size of the lists
Nr = length(r);
Nlam = length(lambda);
Nth = length(theta);

vals = zeros(Nr, Nlam, Nth, 3);

vals(:,:,:,1) = fevalm(F{1}, r, lambda, theta);
vals(:,:,:,2) = fevalm(F{2}, r, lambda, theta);
vals(:,:,:,3) = fevalm(F{3}, r, lambda, theta);
end
