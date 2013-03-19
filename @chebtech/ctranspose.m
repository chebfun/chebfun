function ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE  chebtech objects are not transposable.
    
    error('CHEBFUN:CHEBTECH:ctranspose:notpossible', ...
        'chebtech objects are not transposable.')
    
end