function InputDimensionAdjustment = getInputDimensionAdjustment(L)

if ( isa(L, 'linop') )
    InputDimensionAdjustment = max(getDiffOrder(L), [], 1);
    InputDimensionAdjustment = max(InputDimensionAdjustment, 0);
    InputDimensionAdjustment = repmat(InputDimensionAdjustment, size(L, 1), 1);
else
    InputDimensionAdjustment = 0;
end

end