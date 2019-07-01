function L = laplacian(F)
%LAPLACIAN Vector Laplacian of a BALLFUNV.
%   LAPLACIAN(F) returns a BALLFUNV representing the vector Laplacian of F.
% 
% See also BALLFUN/LAPLACIAN

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( F )
    L = ballfunv();
    return
end

% Compute the vector laplacian
L = [laplacian(F.comp{1}); laplacian(F.comp{2}); laplacian(F.comp{3})];
end
