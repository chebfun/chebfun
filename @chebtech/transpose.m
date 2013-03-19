function transpose(f) %#ok<*INUSD>
%TRANSPOSE  chebtech objects are not transposable.
    
    error('CHEBFUN:CHEBTECH:transpose:notpossible', ...
        'chebtech objects are not transposable.')
    
end