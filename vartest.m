nsamp = 1000;

L = 1;

f = randnfun(.05,[0,L],'norm',nsamp);  
b = cumsum(f);
s = sum(b.^2,2)/nsamp;
plot(s), drawnow
shg

f = randnfun(.05,[0,L],'norm','complex',nsamp);  
b = cumsum(f);
s = sum(conj(b).*b,2)/nsamp;
hold on, plot(s), hold off

legend('real','complex','location','southeast')
