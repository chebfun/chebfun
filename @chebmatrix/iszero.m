function isz = iszero(f)
%ISZERO   Returns a matrix, with entry 1 in the (i, j) positions if the
%         corresponding entry in the cell-array of 1 has .iszero = 1, and 0
%         otherwise.
%
%   This function is used for linearity detection in CHEBOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fBlocks = f.blocks;
isz = zeros(size(fBlocks));

% Loop through all elements.
for j = 1:numel(fBlocks);
    if ( isa(fBlocks{j}, 'linBlock') )
        isz(j) = fBlocks{j}.iszero;
    end
end

end
