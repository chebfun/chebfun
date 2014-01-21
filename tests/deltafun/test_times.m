% Test file for @deltafun/times.m.

function pass = test_times(pref)

% if (nargin < 1)
%     pref = chebpref();
% end
%%
a = -4; b = 4;
pTol = deltafun.pref.deltafun.proximityTol;
dTol = deltafun.pref.deltafun.deltaTol;

f1 = fun.constructor(@(x) exp(sin(x)), [a, b]);
d1 = .5*[0 1 0;
         0 1 0;
         1 0 1; ];
l1 = sort(.9*(a + (b-a)*rand(1,3)));

f2 = fun.constructor(@(x) exp(cos(x)), [a, b]);
d2 = .8*[1 0 1;
         0 1 0;
         0 1 0; ];
l2 = l1 + .0010751808168034;

df1 = deltafun(f1, d1, l1);
df2 = deltafun(f2, [], []);

pass(1) = isempty(deltafun() .* deltafun());
pass(2) = isempty(deltafun() .* df1) && isempty(df1 .* deltafun());

%%
s = df1 .* df2;
s.impulses
pass(3) = iszero(s.funPart - f1.*f2);
pass(4) = norm(s.location - sort(l1), inf) == 0;
c1 = [ feval(diff(f2,2), l1(1));
       -2*feval(diff(f2,1), l1(1));
       feval(diff(f2,0), l1(1));    ];
c2 = [ feval(diff(f,0), l1(2));
       feval(diff(f,1), l1(2));
       0;                       ];
       
c3 = [feval(diff(f2,2), l1(3));
      -2*feval(diff(f2,1), l1(3));
      feval(diff(f2,0), l1(3))];

deltas1 = .5*[c1, c2, c3]

         %%
         delta2 = [ 
pass(5) = norm(s.impulses - A, inf) == 0 && ...
    norm(s.location - sort(l1), inf) == 0 && ...
    iszero( f1+f2 - s.funPart);

l = rand(1, 5);
d = rand(3, 5);
[sl, idx] = sort(l);
A = d(:, idx);
s = deltafun(f1, d, l) + deltafun(f1, A, sl);
pass(4) = norm(s.impulses - 2*A, inf) == 0 && ...
    norm(s.location - sl, inf) == 0 && ...
    iszero( 2*f1 - s.funPart);

df1 = deltafun(f1, [1 2 3 4], [-.25 .5 -.5 -.8]);
df2 = deltafun(f2, [1 2 3 4], [-.25 .6 -.5 -.7]); 
s = df1 + df2;
pass(5) = norm(s.location - [-.8 -.7 -.5 -.25 .5 .6], inf) == 0 && ...
    norm(s.impulses - [4 4 6 2 2 2], inf) == 0;
end