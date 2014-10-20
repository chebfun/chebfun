function pass = test_deal(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

x = chebfun(@(x) x);
f = x;
g = sin(x);
D = operatorBlock.diff;
A = [D, f; g', 2];

a = deal(A);

% Check to make sure single output works. This only checks to make sure it got
% the operator; it doesn't check that the operator has not changed. (Cannot do
% that since there is no `==` for operatorBlocks.)
sizea = size(a);
exact = size(D);
pass(1) = all(sizea == exact);

% Check to make sure multiple outputs work (there is a row-major ordering).
[a,b,c,d] = deal(A);
sizea = size(a);  exacta = size(D);
sizeb = size(b);  exactb = size(f);
sizec = size(c);  exactc = size(g');
sized = size(d);  exactd = size(2);
compare = [ sizea == exacta;
            sizeb == exactb;
            sizec == exactc;
            sized == exactd ];
pass(2) = all(compare(:));

% Make sure the error works.
pass(3) = 0;
try
    [a,b,c,d,e] = deal(A);
catch ME
    if (ME.identifier == 'CHEBFUN:CHEBMATRIX:deal:tooManyOutputs')
        pass(3) = 1;
    end
end


end
