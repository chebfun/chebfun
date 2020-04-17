%% AAA approximation of periodic functions
% Peter J. Baddoo, April 2020

%%
% (Chebfun example approx/AAAtrigApprox.m)
% [Tags: #rational, #AAA, #periodic]
%
% This example introduces the periodic extension of the AAA algorithm, AAAtrig.
%
%% 1. Introduction
%
% The AAA algorithm is the state-of-the-art in rational approximation [1]. 
% It is the most general form of rational approximation available in
% Chebfun, and can be applied to very general domains in the complex
% plane. The success of AAA can broadly be attributed to two
% features of the algorithm: greedy selection of support points, and a
% barycentric representation of the approximant.
% These key features of the AAA can be also be used to approximate periodic
% functions [2].
% For example, whilst the code |aaa| returns a rational approximant in
% barycentric form:
%
% $$\centering r(z) = \sum_{j=1}^m \frac{w_j f_j}{z - z_j} \Bigg/
% \sum_{j=1}^m \frac{w_j}{z - z_j}, $$
%
% the code |aaatrig| returns a $2\pi$-periodic rational approximant in
% _trigonometric_ barycentric form [3]:
%
% $$r(z) = \sum_{j=1}^m {w_j f_j} \textrm{cst}\left(\frac{z - z_j}{2}
% \right) \Bigg/ \sum_{j=1}^m w_j \textrm{cst} \left(\frac{z - z_j}{2}
% \right).$$
%
% Either $\textrm{cst} = \cot$ or $\textrm{cst} = \csc$, 
% depending on the application.
% We refer to these two representations of $\textrm{cst}$ as "even" and "odd" respectively; 
% the default setting is to take an odd approximant.
% 
% In the remainder of this example we will trial the new code |aaatrig| on a
% few test cases.
%
%% 2. Real interval
% As a simple example, let's approximate $\tan(x/2-i)$ sampled uniformly on the interval $[0,2\pi]$:
%
Z = linspace(0,2*pi,1e3);
f = @(t) tan(t/2-1i);
F = f(Z);
[rtrig,pol,res,zer] = aaatrig(F,Z);
disp('        zeros              poles             residues')
disp([zer pol res])
%%
% The code |aaatrig| returns a function handle for the rational approximant, as 
% well as  the poles, zeros, residues, and other
% information about the approximant.
% These poles and zeros are repeated with period $2 \pi$.
% The algorithm correctly identifies the zeros of $f$ at even multiple of $\pi$ and
% the poles at odd multiples of $\pi$, with a shift of $2i$.
%
%% 3. Complex domains
% The domain of approximation need not be the interval $[0,2\pi)$, as in
% other periodic rational interpolants such as |trigratinterp|.
% For example, the domain could be the unit circle, as in the following
% example:
Z = exp(1i*linspace(-pi,pi,1e3)');
f = @(t) exp(1i*tan(t/2));
F = f(Z);
[rtrig, pol, res, zer, zj, fj, wj, errtrig] = aaatrig(F,Z);
[r    , ~, ~, ~, ~, ~, ~, err    ] = aaa(F,Z);
LW = 'LineWidth'; FS = 'FontSize';
semilogy(err,LW,3);
hold on
semilogy(errtrig,LW,3);
legend({'aaa','aaatrig'},FS,20);
grid on; set(gca,'YMinorGrid','off')
ylabel('max error')
hold off

%%
% The error curves show that |aaatrig| performs competitively with |aaa|,
% requiring only a type (7,7) approximation to achieve accuracy of order
% $10^{-14}$.
%
% Once the support points and weights have been calculated, the
% approximant can be differentiated using |diffbarytrig|:
%
dtrig = diffbarytrig(Z, zj, fj,wj);
exact = 1i/2*sec(Z/2).^2.*f(Z); %Exact derivative
disp(['The derivative error is ', num2str(norm(dtrig - exact,inf)), '.']);
plot(angle(Z),real(exact),LW,3);
hold on
plot(angle(Z),real(dtrig),'--',LW,3) % Approximate derivative
legend({'exact derivative','approximate derivative'},FS,12,'Location','southeast');
hold off

%% 4. Unbounded domains
% Rational functions typicall perform better than polynomials when working
% with unbounded domains [4]. Since the previous example we can see
% that the trigonometric approximant remains accurate far away from the
% sample points in the complex plane. For example, evaluating the function
% and its approximants at 500i gives
disp('     function value       aaatrig              aaa')
disp([f(500i), rtrig(500i), r(500i)])

%%
% Looking at the error throughout a period window in the complex plane shows further that
% the trigonometric approximant performs better away from the sample
% points.
x = linspace(-pi,pi); y = linspace(-10,10);
z = x + 1i*y';
subplot(1,2,1);
pcolor(x,y,log(abs(f(z)-rtrig(z)))); hold on
plot(Z,'r.',LW,2); hold off
caxis([-15,0]), title('aaatrig error')
shading interp; axis tight; axis equal
subplot(1,2,2);
pcolor(x,y,log(abs(f(z)-r(z)))); hold on
plot(Z,'r.',LW,2); hold off
caxis([-15,0]), title('aaa error')
shading interp; axis tight; axis equal
cb=colorbar;cb.Location = 'eastoutside'; 
cb.Label.String = 'log of error';

%%
% Moreover, a new feature of |aaatrig| is that values at $\pm i \times
% \infty$ can also be specified. For example,
Z = linspace(0,2*pi,2e3);
F = gamma(2+sin(Z));
F = [F, -(1+sqrt(5))/2,     pi];
Z = [Z,         1i*Inf, -1i*Inf];
rtrig  = aaatrig(F,Z);

%%
% Now the approximant takes the desired values at infinity.
disp(rtrig([1i*Inf,-1i*Inf]))

%% References 
% [1] Y. Nakatsukasa, O. Sète, and L. N. Trefethen, 
% “The AAA Algorithm for Rational Approximation,” 
%  _SIAM J. Sci. Comput._, vol. 40, no. 3, pp. A1494–A1522, 2018.
%
% [2] P. J. Baddoo,
% “AAA rational approximation of periodic functions,” 
%  _In Prepration_, 2020.
%
% [3] J.-P. Berrut, 
% “Baryzentrische Formeln zur Trigonometrischen Interpolation (I),” 
% _ZAMP Zeitschrift für Angew. Math. und Phys._, vol. 35, no. 1, pp. 91–105, 1984.
%
% [4] L. N. Trefethen, _Approximation Theory and Approximation
% Practice, extended edition_, SIAM, 2020.

