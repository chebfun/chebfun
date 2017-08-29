function coreplot(f, varargin)
%COREPLOT   A scatter3 plot of the core tensor of a CHEBFUN3 that 
%   visualizes the logarithm of the magnitude of each entry.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    return
end

T = f.core;
len = size(T);
if ( ndims(T) < 3 )
    % This happens, e.g. for chebfun3 of a univariate function
   len = [len 1];
end
    
[ind1, ind2, ind3] = ndgrid(1:len(1), 1:len(2), 1:len(3));
scatter3(ind1(:), ind2(:), ind3(:), [], log10(abs(T(:))))
colorbar

title('logarithm of the magnitude of entries in the core tensor')

end