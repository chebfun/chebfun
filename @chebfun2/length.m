function [out1, out2] = length(F)

if ( nargout == 1 )
    out1 = length(F.pivotValues);
else
    out1 = length(F.rows);
    out2 = length(F.cols);
end

end