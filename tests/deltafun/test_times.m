% Test file for @deltafun/times.m.

function pass = test_times(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
% Get the tolerance:
dTol = pref.deltaPrefs.deltaTol;

d = deltafun(bndfun(@sin), struct('deltaMag', 1, 'deltaLoc', 0));
pass(1) = isempty(deltafun() .* deltafun());
pass(2) = isempty(deltafun() .* d) && isempty(d .* deltafun());

f = fun.constructor(@(x) exp(-x));
g = bndfun(@sin);
df1 = deltafun(f, []);
df2 = deltafun(g, struct('deltaMag', [0; 0; 0; 0; 1], 'deltaLoc', 0));
s = df1.*df2;
pass(3) = norm(s.deltaMag - [1, 4, 6, 4, 1].', inf) < dTol;

f = fun.constructor(@(x) exp(x));
g = bndfun(@sin);
df1 = deltafun(f, []);
df2 = deltafun(g, struct('deltaMag', [0; 0; 0; 1], 'deltaLoc', 0));
s = df1.*df2;
pass(4) = norm(s.deltaMag - [-1, 3, -3, 1].', inf) < dTol;

a = -4;
b = 4;

f1 = fun.constructor(@(x) exp(sin(x)), struct('domain', [a, b]));
d1 = .5*[0 1 0;
         0 1 0;
         1 0 1; ];
l1 = sort(.9*(a + (b-a)/2*rand(1,3)));

f2 = fun.constructor(@(x) exp(cos(x)), struct('domain', [a, b]));

df1 = deltafun(f1, struct('deltaMag', d1, 'deltaLoc', l1));
df2 = deltafun(f2, []);

s = df1 .* df2;
pass(5) = iszero(s.funPart - f1.*f2);
pass(6) = norm(s.deltaLoc - sort(l1), inf) == 0;

c1 = [ feval(diff(f2,2), l1(1));
       -2*feval(diff(f2,1), l1(1));
       feval(diff(f2,0), l1(1));    ];

c2 = [ feval(diff(f2,0)-diff(f2,1), l1(2));
       feval(diff(f2,0), l1(2));
       0;                       ];
       
c3 = [feval(diff(f2,2), l1(3));
      -2*feval(diff(f2,1), l1(3));
      feval(diff(f2,0), l1(3))];

deltas1 = .5*[c1, c2, c3];

err = s.deltaMag - deltas1;
pass(7) = norm(err(:), inf) < dTol;
    

d2 = .8*[1 0 1;
         0 1 0;
         0 1 0; ];
l2 =  .0010751808168034 + sort(.9*(b-a)/2*rand(1,3));
df2 = deltafun(f2, struct('deltaMag', d2, 'deltaLoc',l2));
%%
s = df1 .* df2;
pass(8) = iszero(s.funPart - f1.*f2);
pass(9) = norm(s.deltaLoc - sort(union(l1,l2)), inf) == 0;
%%
c1 = [ feval(diff(f1,0), l2(1));
       0;
       0; ];
c2 = [ feval(diff(f1, 2)-diff(f1, 1), l2(2));
       feval(diff(f1, 0)-2*diff(f1, 1), l2(2));
       feval(diff(f1, 0), l2(2));  ];
       
c3 = [feval(diff(f1,0), l2(3));
      0;
      0; ];

deltas2 = .8*[c1, c2, c3];

err = [deltas1, deltas2] - s.deltaMag;
pass(10) = norm(err(:), inf) < dTol;
%%
end
