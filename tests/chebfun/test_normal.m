function pass = test_normal( pref )
% Test normal vector to complex chebfuns.

if ( nargin == 0 )
    pref = chebfunpref;
end

tol = pref.chebfuneps;

%% #2268
g = chebfun(@(t) 0.3*cos(t) + 1i*sin(t), [0, 2*pi]);
ntrue = -1i*diff(g);
ntrue = ntrue ./ abs(ntrue);
ntrue = [real(ntrue), imag(ntrue)];
n = normal(g, 'unit');

pass(1) = ( norm(n - ntrue) < tol );

end
