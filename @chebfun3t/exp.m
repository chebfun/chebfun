function f = exp(f) 
%EXP  Exponential of a CHEBFUN3T
%   EXP(F) returns the exponential of a CHEBFUN3T. 

% Empty check:
if ( isempty( f ) ) 
    return 
end 

op = @(x,y,z) exp( feval(f, x, y, z ) ); % Resample.
f = chebfun3t( op, f.domain );       % Call constructor.

end
