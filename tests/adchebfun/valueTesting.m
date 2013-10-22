function error = valueTesting(f)

% Seed random generator to ensure same values.
seedRNG(6179);

% Generate an arbitrary CHEBFUN to evaluate the function at
u = chebfun(rand(15,1)) + 1;

% Construct a corresponding ADCHEBFUN
v = adchebfun(u);

% Evaluate f at both u and v. FU will be a CHEBFUN, while FV will be an
% ADCHEBFUN
fu = f(u);
fv = f(v);

% We should expect that the FUNC field of FV matches the FU
error = norm(fu-fv.func,'inf');

end