function F = fliplr(F)
%FLIPLR   Flip the columns of a CHEBMATRIX.

F.blocks = fliplr(F.blocks);

end