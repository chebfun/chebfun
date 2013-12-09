function [out1, out2] = length(F)
% LENGTH   Rank of the chebfun2

if ( nargout <= 1 )
    out1 = length(F.pivotValues);
else
    out1 = length(F.rows);
    out2 = length(F.cols);
end

end