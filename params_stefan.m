function p = params_stefan

p.eps = sqrt(pi*((3.5e-6)^2 + (2e-5)^2));
%p.eps = 1e-3;

p.length = 0.5;
p.r2 = 0.45*p.eps;

p.cw = 4184;
p.Dw = 0.58/1000;

p.Tc = 273.15;
p.TempOut = p.Tc + 10;
p.TimeMax = 5e6;


p.kw = 0.58;
p.lambda = 3.34e5;
p.rhow = 1000;


%Initial values
p.Rv = 0.5557439*p.eps;    
p.Rf = 0.0972552*p.eps;
p.TempInit = p.Tc;
p.sig0 = p.Rf/sqrt(2);
p.siw0 = p.Rf;
p.r0 = 0.3*p.Rv;
p.U0 = 0;


%% Diffusion tensor

%Values for homogenized diffusion coeffizient not considering gas
r=p.eps*[0 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.47 0.49 0.5];
Pr=[1 0.99937 0.98442 0.93909 0.86796 0.77672 0.67163 0.55849 0.441496 0.32209 0.19649 0.14064 0.07274 0];
p.function_PT = @(x) p.Dw*spline(r,Pr,x);
%Values for homogenized diffusion coeffizient considering gas
sig_difftensor=p.eps*[-1000 0.1 0.2 0.3 0.369 1000];
r_difftensor=p.eps*[-1000 0 0.05 0.1 0.15 0.2 1000];
[SIG_difftensor,R_difftensor]=meshgrid(sig_difftensor,r_difftensor);
Values_difftensor = [1.3147e-7   1.3147e-7   1.0874e-7   7.8188e-8   5.5503e-8   5.5503e-8;...
                    1.3147e-7    1.3147e-7   1.0874e-7   7.8188e-8   5.5503e-8   5.5503e-8;...
                    2.7538e-7    2.7538e-7   2.4209e-7   1.9084e-7   1.5264e-7   1.5264e-7;...
                    7.0615e-7    7.0615e-7   6.3465e-7   5.3134e-7   4.5363e-7   4.5363e-7;...
                    1.4246e-6    1.4246e-6   1.2966e-6   1.1059e-6   9.0573e-7   9.0573e-7;...
                    2.4200e-6    2.4200e-6   2.2084e-6   0           0           0;...
                    2.4200e-6    2.4200e-6   2.2084e-6   0           0           0];
 Values_difftensor_melted = [7.3925e-7   7.3925e-7   2.5050e-6   5.4733e-6   8.2459e-6   8.2459e-6;...
                             7.3925e-7   7.3925e-7   2.5050e-6   5.4733e-6   8.2459e-6   8.2459e-6;...
                             8.8625e-7   8.8625e-7   2.6534e-6   5.6217e-6   8.3943e-6   8.3943e-6;...
                             1.3423e-6   1.3423e-6   3.1047e-6   6.0730e-6   8.6585e-6   8.6585e-6;...
                             2.0741e-6   2.0741e-6   3.8399e-6   6.7266e-6   8.9126e-6   8.9126e-6;...
                             3.0965e-6   3.0965e-6   4.8622e-6   0           0           0;...
                             3.0965e-6   3.0965e-6   4.8622e-6   0           0           0];
 p.function_difftensor = @(x,y) p.cw*interp2(SIG_difftensor,R_difftensor,Values_difftensor,x,y);
 p.function_difftensor_melted =  @(x,y) p.cw*interp2(SIG_difftensor,R_difftensor,Values_difftensor_melted,x,y);

%% Preparations Global

p.LG = p.length/2;
p.NG = 50;
p.hG = p.LG/(p.NG+1);
h12 = p.LG/(2*p.NG+1);
p.langG = h12:2*h12:p.LG-p.hG;
A1 = -diag(2*ones(p.NG,1)) + diag(ones(p.NG-1,1),1) + diag(ones(p.NG-1,1),-1);
A1(1,1) = -1;
A2 = diag(1./p.langG(1:p.NG-1),1) - diag(1./p.langG(2:p.NG),-1);
A2(1,1) = -1/p.langG(1);
p.AG = A1/p.hG^2 + A2/2/p.hG;
p.bG = [zeros(p.NG-1,1);p.TempOut/p.hG*(1/p.hG+1/2/p.langG(p.NG))];

%% Preparations Local

p.LL = 1;
p.NL = 4;  % only the inner points;
p.hL = p.LL/(p.NL+1);
p.langL = p.hL*(0:p.NL+1);
p.A3 = -diag(2*ones(p.NL,1)) + diag(ones(p.NL-1,1),1) + diag(ones(p.NL-1,1),-1);


%% Cell parameters
p.Lf = 1e-3;        % m
p.Lv = 0.5e-3;      % m
p.A = 2*pi*p.Rf*p.Lf;  % m^2
p.Vf = pi*p.Rf^2*p.Lf; % m^3
p.Vv = pi*p.Rv^2*p.Lv;  % m^3

p.cs = 58.4;        % mol/m^3
p.g = 9.81;         % m/s^2
p.H = 0.0274;       % - 
p.K = 1.98e-14;     % m/s
p.Mg = 0.029;       % kg/mol
p.Nf = 16;           % -
p.R = 8.314;        % J/molK
p.W = 0.101145395*p.eps; % m
p.rhoi = 917;       % kg/m^3
p.sigmaw = 76e-3;   % N/m
p.sigmas = p.sigmaw;  %75.6e-3; % N/m

%Initial values
p.pfg0 = 2e5;   % N/m^2
p.pvg0 = 1e5;    % N/m^2

p.pfw0 = p.pfg0 - p.sigmaw/p.sig0;  
p.Vvg0 = pi*p.Lf*p.r0^2;
p.rhovg0 = p.pvg0*p.Mg/p.R/p.Tc;
p.pvw0 = p.pvg0 - p.sigmas/p.r0;
end
