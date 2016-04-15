function f = mrdivide(f, g)
%/   Right scalar divide for CHEBFUN3 objects.
%
%    F/C divides the CHEBFUN3 object F by a scalar C.

if ( isa(g, 'double') )
    f.core = f.core / g;
else
    error('CHEBFUN:CHEBFUN3:mrdivide:mrdivide', ...
        'Not supported. Did you mean ./ ?');
end

end