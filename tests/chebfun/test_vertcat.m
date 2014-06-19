function pass = test_vertcat(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% [CHEBFUN ; DOUBLE]
x = chebfun('x', pref);
f = [x ; 1];
pass(1) = isa(f, 'chebmatrix') && all(size(f.blocks)==[2,1]) && ...
    isnumeric(f.blocks{2,1});

% [CHEBFUN ; CHEBFUN]
f = [x ; x];
pass(2) = isa(f, 'chebmatrix') && all(size(f.blocks)==[2,1]) && ...
    isa(f.blocks{2,1}, 'chebfun');

% Cast double to CHEBFUN
f = [x 1 ; x 1];
pass(3) = isa(f, 'chebmatrix') && all(size(f.blocks)==[2,2]) && ...
    isa(f.blocks{2,1}, 'chebfun') && isa(f.blocks{1,2}, 'chebfun');

% Array-valued CHEBFUN
f = [x x ; x x];
pass(4) = isa(f, 'chebmatrix') && all(size(f.blocks)==[2,2]) && ...
    isa(f.blocks{2,1}, 'chebfun') && isa(f.blocks{1,2}, 'chebfun');

% Quasimatrix
f = [x abs(x) ; x x];
pass(4) = isa(f, 'chebmatrix') && all(size(f.blocks)==[2,2]) && ...
    isa(f.blocks{2,1}, 'chebfun') && isa(f.blocks{1,2}, 'chebfun');

% Row CHEBFUNs:
x = x';
f = [x ; x];
pass(5) = isa(f, 'chebfun') && numColumns(f) == 2;

% Incorrect transpose state:
try 
    f = [x ; x.'];
    pass(6) = false;
catch
    pass(6) = true;
end

end