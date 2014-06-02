function pass = test_flip(pref)
% GBW, 7 May 2014

%% 
% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

d = [-2 2];
x = chebfun('x',d);
a11 = sign(x);
a12 = cos(pi*x)+1i*sin(x);
a13 = operatorBlock.mult(x.^2);
a21 = pi;
a22 = 2*pi;
a23 = functionalBlock.sum(d);
a31 = atan(x);
a32 = abs(a11);
a33 = operatorBlock.diff(d);

A = [a11,a12,a13;a21,a22,a23;a31,a32,a33];
B = flipud(A);

count = 1;
for j=1:3
    pass(count) = isequal(A{1,j},B{3,j});
    count = count + 1;
end

for j=1:3
    pass(count) = isequal(A{2,j},B{2,j});
    count = count + 1;
end

B = fliplr(A);

for j=1:3
    pass(count) = isequal(A{j,1},B{j,3});
    count = count + 1;
end

for j=1:3
    pass(count) = isequal(A{j,2},B{j,2});
    count = count + 1;
end

end
