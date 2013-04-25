function ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE  FUN objects are not transposable.
    
error('CHEBFUN:FUN:ctranspose:notpossible', ...
    'FUN objects are not transposable.')
    
end
