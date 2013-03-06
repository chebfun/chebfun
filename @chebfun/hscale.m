function h = hscale(f)

h = norm(f.domain, inf);
if ( isinf(h) )
    h = 1;
end