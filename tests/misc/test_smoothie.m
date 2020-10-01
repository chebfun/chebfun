function pass = test_smoothie(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

f = smoothie;
s = std(f);
pass(1) = (s>0.2) & (s<5);

f = smoothie('trig');
s = std(f);
pass(2) = (s>0.2) & (s<5);

f = smoothie([-.1 .1]);
s = std(f);
pass(3) = (s>0.1) & (s<5);

f = smoothie([-3 3]);
s = std(f);
pass(4) = (s>0.2) & (s<5);

f = smoothie('complex');
s = std(f);
pass(5) = (s>0.2) & (s<5);

f = smoothie('trig','complex',[-7,-6]);
s = std(f);
pass(6) = (s>0.2) & (s<5);

f = smoothie(3);
pass(7) = (size(f,2) == 3);

end
