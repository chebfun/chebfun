%% Fourier Transforms using Contour Integrals
% Mohsin Javed, 3rd July 2013.
 
%%
% (Chebfun example complex/FourierTransforms.m)
% [Tags: #fouriertransform, #contourintegral, #pole, #residue]
LW = 'linewidth'; lw = 1.5; 
MS = 'markersize';
%% Fourier Transform of a Rational Function
% The Fourier Transform of a function $f$ is defined as
% $$ F(w) = \frac{1}{2\pi}\int_{-\infty}^{\infty}f(z)e^{-iwz} dz. $$
%%
% Let us consider the rational function
% $$ f(z) = \frac{2a}{a^2+z^2}. $$
% Then, for $w > 0$, we form a closed contour running 
% along the real axis from $-R$ to $R$ and then 
% along a semi-circle $\Gamma_R$ of radius $R$ in the
% lower half plane. This contour encloses the pole
% of $f$ at $z=-ia$. Also, the contour is oriented negatively,
% therefore, by the residue theorem [1]:
% $$ \int_{-R}^{R}f(z)e^{-iwz} dz + \int_{\Gamma_R}f(z)e^{-iwz} dz = (-2\pi i) Res(f(z)e^{-iwz}, -ia). $$
R = 5;
x = chebfun('x', [-R,R]) + eps*1i;
z = chebfun(@(t) R*exp(-1i*t), [0, pi] );
plot(x), hold on
plot(z), axis equal
xlim( [-R-1, R+1] )
%%
% Note that for $w < 0$, we can use a similar 
% contour in the upper half plane.
%%
% It can be shown (see, for example [1]) that 
% as $R \to \infty$, the integral along 
% $\Gamma_R$ vanishes and hence the value 
% of the Fourier transform at $w$ is given as 
% $$ F(w) = \frac{1}{2 \pi}(-2\pi i) Res(f(z)e^{-iwz}, -ia). $$
% The residue above can be computed
% very easily and accurately in Chebfun.
% We first form a circular contour
% lying in the lower half plane 
% and enclosing only the pole at $z=-ia$.
a = 1.5;
r = a/2;
z = chebfun(@(t) -a*1i+ r*exp(1i*t), [0, 2*pi] );
f = 2*a./(a^2+z.^2);
plot(z), hold on
plot( 0, -a, 'xr', MS, 12 ), hold off
%%
% We then define a function handle $\verb|Fh|$ 
% for the value of the Fourier transform 
% at a point $w$. 
Fh = @(w) -1/(2*pi)*sum(f.*exp(-1i*w.*z).*diff(z));
%%
% Finally, we form a chebfun for 
% the Fourier transform of 
% $f$.
F = chebfun(Fh, [0, 10], 'vectorize')

%%
% Ideally $F(w)$ should be real.
norm(imag(F), inf)
%%
% Let's have a look 
% at the Fourier 
% transform of $f$.
F = real(F);
plot(F, LW, lw), hold on
title( 'The Fourier Transform of 2a/(x^2+a^2)' )
%% 
% The Fourier transform for this particular $f$ is 
% known in closed form:
% $$ F(w) = e^{-a|w|}. $$
Fexact = chebfun( @(w) exp(-a*w), [0, 10 ] )
plot( Fexact, '.r', MS, 15), hold off
norm(F-Fexact, inf )

%% Fourier Transform of $1/x$.
% We saw in the previous example that 
% $|f|$ decayed like $O(1/R^2)$ on $\Gamma_R$, which
% guaranteed the vanishing of the contour integral
% on $\Gamma_R$ as $R \to \infty$. However, a result
% known as Jordan's lemma shows
% that even if the function merely goes 
% to zero uniformly as $R \to \infty$, 
% the contour integral
% on $\Gamma_{R}$ still vanishes and we can apply
% the same technique [1].

%% 
% The Fourier transform of $1/x$ is defined as a principal
% value integral:
% $$ F(w) = \frac{1}{2\pi}PV\int_{-\infty}^{\infty}\frac{1}{z}e^{-iwz} dz. $$
% Repeating the same procedure, we first 
% form a circular contour 
% around the origin.
R = 0.2;
z = chebfun(@(t) R*exp(1i*t), [0, 2*pi] );
f = 1./z;
%%
% We then form a function handle for 
% the value of the Fourier transform 
% of $1/z$ at a point $w$. Since the pole
% now lies on the original contour, following
% standard techniques of principal value integrals,
% we only pick half of the residue at the origin.
Fh = @(w) (1/2)*(-1/(2*pi))*sum(f.*exp(-1i*w.*z).*diff(z));
%%
% Finally, we form a chebfun for 
% the Fourier transform of 
% $f$.
F = chebfun(Fh, [0, 4], 'vectorize')

%%
% Since $f$ is odd, $F(w)$ should be purely imaginary.
% Let's have a look 
% at the Fourier 
% transform of $f$.
norm(real(F), inf)
imagF = imag(F)
plot(imagF, LW, lw)
%% 
% The Fourier transform for this particular $f$ is 
% known in closed form:
% $$ F(w) = -\frac{i}{2}\mbox{sign}(w). $$
norm(F-(-1i/2), inf )
%%
% For this function, we can utilize Chebfun's 
% splitting on mode to
% construct the Fourier transform for
% positive and negative values of $w$ simultaneously. For $w < 0$, 
% the contour $\Gamma_R$ should be closed in the upper
% half plane in order to apply Jordan's lemma. This results in 
% a flipping of the sign for 
% the residue at the origin when $w<0$. The Chebfun 
% representation of $F(w)$ can now be obtained as follows.
Fh = @(w) -1/(2*pi)*(1/2)*sum(f.*exp(-1i*w.*z).*diff(z)).*sign(w);
F = chebfun(Fh, [-1, 1], 'vectorize', 'splitting', 'on')
%%
% Again, the Fourier transform should ideally be
% purely imaginary but we can expect rounding errors. 
% Let's check the accuracy
% of our computed Fourier transform.
norm(real(F), inf )
Fexact = chebfun(@(w) -1i/2*sign(w), [-1,1], 'splitting', 'on' );
norm(F-Fexact, inf )
%%
imagF = imag(F)
plot(imagF, LW, lw), ylim( [-.6, .6] )
title( 'The imaginary part of the Fourier transform of 1/x.' )
%% Fourier Transform of $1/(x-a)$
% We also know from the properties of 
% Fourier transform that if $f(z)$ has
% the Fourier transform $F(w)$, then
% given a real number $a$, the transform 
% of $f(z-a)$ is $e^{-iaw}F(w)$.

%%
% We now construct the Fourier transform of
% $1/(x-a)$ for $a = 2\pi$.
a = 2*pi;
R = 0.2;
z = chebfun(@(t) a + R*exp(1i*t), [0, 2*pi] );
f = 1./(z-a);
Fh = @(w) -1/(2*pi)*(1/2)*sum(f.*exp(-1i*w.*z).*diff(z)).*sign(w);
F = chebfun(Fh, [-1, 1], 'vectorize', 'splitting', 'on')
%%
% We then check the accuracy
% of our computation.
Fexact = chebfun(@(w) -1i/2*sign(w).*exp(-1i*a*w), [-1,1], 'splitting', 'on' );
norm(F-Fexact, inf )
%%
% Here is a plot of the real 
% and imaginary part of $F(w)$.
subplot(1,2,1)
plot(real(F), LW, lw), ylim( [-.6, .6] )
title( 'The real part of the Fourier transform of 1/(x-a).' )
subplot(1,2,2)
plot(imag(F), LW, lw), ylim( [-.6, .6] )
title( 'The imaginary part of the Fourier transform of 1/(x-a).' )

%% References
% 
% [1] Saff, E. B. and Snider, A. D., Fundamentals of Complex Analysis,  
% Prentice-Hall, 2003.