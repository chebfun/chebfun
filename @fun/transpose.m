function transpose(f) %#ok<*INUSD>
%TRANSPOSE   BNDFUN objects are not transposable.

error('CHEBFUN:FUN:transpose:notpossible', ...
    'FUN objects are not transposable.')

end
