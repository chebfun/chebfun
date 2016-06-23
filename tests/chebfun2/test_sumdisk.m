function pass = test_sumdisk( pref ) 
% Test for sumdisk command of a chebfun2 object. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 100*eps;  % intentionally scaled to macheps, not chebfun2eps
j = 1; 
pass = [];

% Load set of test functions

ListLength = 18;

TestFn_Cell = cell(ListLength,5);           
% 1st column function name, 2nd column function itself
% 3rd column is a non-defult domain
% 4th column excat intergal on default domain given by integral 2
% 5th column excat intergal on non-default domain given by integral 2

TestFn_Cell{1,1} = 'Exp fn';
TestFn_Cell{1,2} = @(x,y) exp(-30*(x.^2+y.^2));
TestFn_Cell{1,3} = [0 3 0 3];
TestFn_Cell{1,4} = 0.10471975511965049555;
TestFn_Cell{1,5} = 6.5047850792584266207e-08;

TestFn_Cell{2,1} = 'Bessel fn';
TestFn_Cell{2,2} = @(x,y) besselj(0,10*sqrt(x.^2+y.^2));
TestFn_Cell{2,3} = [-2 0 -2 0];
TestFn_Cell{2,4} = 0.027314731999093719295;
TestFn_Cell{2,5} = 0.0041123909829739284383;

TestFn_Cell{3,1} = 'Cos fn 1';
TestFn_Cell{3,2} = @(x,y) .3*cos(2*(x.^2+y.^2)/sqrt(2));
TestFn_Cell{3,3} = [-2 2 -2 2];
TestFn_Cell{3,4} = 0.65827927025174792774;
TestFn_Cell{3,5} = -0.39064683099739677674;

TestFn_Cell{4,1} = 'Runge fn';
TestFn_Cell{4,2} = @(x,y) 1./(0.1+4*(x.^2+y.^2));
TestFn_Cell{4,3} = [1 2 3 4];
TestFn_Cell{4,4} = 2.916632680833628477;
TestFn_Cell{4,5} = 0.013635312695040419442;

TestFn_Cell{5,1} = 'Cos fn 2';
TestFn_Cell{5,2} = @(x,y) (1.25+cos(5*y))./(6+6*(3*x-1).^2);
TestFn_Cell{5,3} = [-1 2 -2 1];
TestFn_Cell{5,4} = 0.25256231557730668413;
TestFn_Cell{5,5} = 0.49034148647026926104;

TestFn_Cell{6,1} = 'test fn 1';
TestFn_Cell{6,2} = @(x,y) (1/3)*cos(x.*2 - y.^2);   
TestFn_Cell{6,3} = [1 2 1 2];
TestFn_Cell{6,4} = 0.55528384075791314967;
TestFn_Cell{6,5} = 0.128338166277416732;

TestFn_Cell{7,1} = 'test fn 2';
TestFn_Cell{7,2} = @(x,y) (1/3)*sin(x - y) .* cos(x + y);
TestFn_Cell{7,3} = [1 3 0.5 1];
TestFn_Cell{7,4} = -1.5039590903762663915e-15;
TestFn_Cell{7,5} = -0.18366707024222650446;

TestFn_Cell{8,1} = 'test fn 3';
TestFn_Cell{8,2} = @(x,y) x.^2 + y - 0.5;
TestFn_Cell{8,3} = [-0.5 0 -3 0];
TestFn_Cell{8,4} = -0.785398163397448279;
TestFn_Cell{8,5} = -2.2641556429192055688;

TestFn_Cell{9,1} = 'test fn 4';
TestFn_Cell{9,2} = @(x,y) 2*cos(10*x).*sin(10*y) + sin(10*x.*y);
TestFn_Cell{9,3} = [2 4 2 4];
TestFn_Cell{9,4} = -7.4348211066127805474e-17;
TestFn_Cell{9,5} = -0.037050942948795864695;

TestFn_Cell{10,1} = 'even fn';
TestFn_Cell{10,2} = @(x,y) exp(-50*x.^2) + 0.75*exp(-50*y.^2) +0.75*exp(-50*x.^2).*exp(-50*y.^2);
TestFn_Cell{10,3} = [-2 4 -2 4];
TestFn_Cell{10,4} = 0.92002342600094488834;
TestFn_Cell{10,5} = 2.5268118478733239129;

TestFn_Cell{11,1} = 'odd fn';
TestFn_Cell{11,2} = @(x,y) (sin(x) + sin(y)).*(exp(x.^2));
TestFn_Cell{11,3} = [0 1 -1 5];
TestFn_Cell{11,4} = -5.4916678097886341414e-15;
TestFn_Cell{11,5} = 4.9596747009976374088;

TestFn_Cell{12,1} = 'const fn';
TestFn_Cell{12,2} = @(x,y) 1;
TestFn_Cell{12,3} = [-3 3 -1 4];
TestFn_Cell{12,4} = 3.1415926535897922278;
TestFn_Cell{12,5} = 23.561944901923439488;

TestFn_Cell{13,1} = 'linear fn';
TestFn_Cell{13,2} = @(x,y) x + y;
TestFn_Cell{13,3} = [-3 0 0 3];
TestFn_Cell{13,4} = -5.7041142608718533452e-15;
TestFn_Cell{13,5} = -1.9505843168704083829e-14;

TestFn_Cell{14,1} = 'quadratic fn 3';
TestFn_Cell{14,2} = @(x,y) x.^2 + y.^2;
TestFn_Cell{14,3} = [-0.5 0 -3 0];
TestFn_Cell{14,4} = 1.5707963267948974462;
TestFn_Cell{14,5} = 3.4054373491061209478;

TestFn_Cell{15,1} = 'trigfun fn 1';
TestFn_Cell{15,2} = @(x,y) cos(pi*x) + sin(pi*y);
TestFn_Cell{15,3} = [1 3 1 3];
TestFn_Cell{15,4} = 0.56923068635950269112;
TestFn_Cell{15,5} = 0.56923068635950269112;

TestFn_Cell{16,1} = 'trigfun fn 2';
TestFn_Cell{16,2} = @(x,y) cos(pi*x).*sin(pi*y);
TestFn_Cell{16,3} = [0 6 0 4];
TestFn_Cell{16,4} = -2.9264458932044042777e-17;
TestFn_Cell{16,5} = -7.3509610067165625678e-16;

TestFn_Cell{17,1} = 'trigfun fn 3';
TestFn_Cell{17,2} = @(x,y) (cos(2*pi*x)).^2 + sin(pi*y);
TestFn_Cell{17,3} = [1.5 3.5 -1 3];
TestFn_Cell{17,4} = 1.5321636228988491091;
TestFn_Cell{17,5} = 3.0643272457977017709;

TestFn_Cell{18,1} = 'trigfun fn 4';
TestFn_Cell{18,2} = @(x,y) 1./(sin(2*pi*x) + 2);
TestFn_Cell{18,3} = [10 12 -2 2];
TestFn_Cell{18,4} = 1.8200461512756580529;
TestFn_Cell{18,5} = 3.6400923025513156617;

% Test sumdisk using chebfun2 test functions defined on the unit disk
Int_Err_Vec = zeros(36,1);

for I = 1:14

    f = TestFn_Cell{I,2};
   
    Integral_Disc_integral2 = TestFn_Cell{I,4};
    
    f_Chebfun2 = chebfun2(f, [-1 1 -1 1]);
    Integral_Disc_sumdisk = sumdisk(f_Chebfun2);
    
    Int_Err_Vec(j) = Integral_Disc_integral2 - Integral_Disc_sumdisk;
    pass(j) = abs(Int_Err_Vec(j) < tol*vscale(f_Chebfun2));
    j = j+1;
end

% Test sumdisk using trigfun2 test functions defined on the unit disk

Int_Err_Vec = zeros(14,1);

for I = 15:18
    f = TestFn_Cell{I,2};
  
    Integral_Disc_integral2 = TestFn_Cell{I,4};
    
    f_Chebfun2 = chebfun2(f, 'trig');
    Integral_Disc_sumdisk = sumdisk(f_Chebfun2);
    
    Int_Err_Vec(j) = Integral_Disc_integral2 - Integral_Disc_sumdisk;
    pass(j) = abs(Int_Err_Vec(j) < tol*vscale(f_Chebfun2));
    j = j+1;
end

% Test sumdisk using chebfun2 test functions defined on the non-default
% domain
Int_Err_Vec = zeros(14,1);

for I = 1:14

    fdomain = TestFn_Cell{I,3};
    f = TestFn_Cell{I,2};
    
    Integral_Disc_integral2 = TestFn_Cell{I,5};
    
    f_Chebfun2 = chebfun2(f, fdomain);
    Integral_Disc_sumdisk = sumdisk(f_Chebfun2);
    
    Int_Err_Vec(j) = Integral_Disc_integral2 - Integral_Disc_sumdisk;
    pass(j) = abs(Int_Err_Vec(j) < tol*vscale(f_Chebfun2));
    j = j+1;
end
 
% Test sumdisk using trigfun2 test functions defined on the non-default
% domain

Int_Err_Vec = zeros(14,1);

for I = 15:18

    fdomain = TestFn_Cell{I,3};
    f = TestFn_Cell{I,2};    
    
    Integral_Disc_integral2 = TestFn_Cell{I,5};
    
    f_Chebfun2 = chebfun2(f, fdomain, 'trig');
    Integral_Disc_sumdisk = sumdisk(f_Chebfun2);
    
    Int_Err_Vec(j) = Integral_Disc_integral2 - Integral_Disc_sumdisk;
    pass(j) = abs(Int_Err_Vec(j) < tol*vscale(f_Chebfun2));
    j = j+1;
end

Max_Err = max(abs(Int_Err_Vec));

end