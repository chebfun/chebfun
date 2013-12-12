function B = mpower(A, pow)

if ( pow ~= round(pow) || pow < 0 )
    error('Power must be a positive integer.')
end

% Create an "identity" chebmatrix for the given variable types.
B = blockEye(A);
       
for i = 1:pow
    B = B*A;
end

end
