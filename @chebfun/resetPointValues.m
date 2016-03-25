function f = clearPointValues(f)
%   F = CLEARPOINTVALS(F) sets the .POINTVALS(j) to the average of the right and
%   left limits of its neighbouring funs for interior breaks and the limits from
%   the left and right for the POINTVALS(1) and VALS(end), respectively.
for k = 1:numel(f)
    f(k).pointValues = chebfun.getValuesAtBreakpoints(f(k).funs);
end

end


