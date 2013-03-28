function transpose(f) %#ok<*INUSD>
%TRANSPOSE   CHEBTECH objects are not transposable.
    
    error('CHEBFUN:CHEBTECH:transpose:notpossible', ...
        'CHEBTECH objects are not transposable.')
    
end
