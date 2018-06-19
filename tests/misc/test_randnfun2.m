function pass = test_randnfun2( pref )

if ( nargin == 0 ) 
    pref = chebfunpref();
end

rng(0)
f = randnfun2(.1);
pass(1) = (abs(mean2(f.^2)-1) < .1);
pass(2) = (abs(mean2(f)) < .1);

rng(0)
g = randnfun2(.1,[2 4 -1 1]);
pass(3) = abs( f(.5,.5) - g(3.5,.5) ) < 1e-12;

rng(0) 
h = randnfun2(.2,[-2 2 0 4]);
pass(4) = abs( f(.5,.5) - h(1,3) ) < 1e-12;

f = randnfun2(1e6);
pass(5) = norm(diff(f),'fro') < 1e-4;

rng(0)
f = randnfun2(.1,'trig');
pass(6) = (abs(mean2(f.^2)-1) < .1);
pass(7) = (abs(mean2(f)) < .1);

rng(0)
g = randnfun2(.1,'trig',[2 4 -1 1]);
pass(8) = abs( f(.5,.5) - g(3.5,.5) ) < 1e-12;

rng(0) 
h = randnfun2('trig',.2,[-2 2 0 4]);
pass(9) = abs( f(.5,.5) - h(1,3) ) < 1e-12;

f = randnfun2(6,'trig');
pass(10) = norm(diff(f),'fro') == 0;

% test that 'big', as well as 'norm', work 
rng(0), f1 = sqrt(10)*randnfun2(.1);
rng(0), f2 = randnfun2(.1,'big');
rng(0), f3 = randnfun2(.1,'norm');
pass(11) = norm(f1-f2) == 0;
pass(12) = norm(f2-f3) == 0;



