% Run this file to newly calculate the solution of the thawing process. The solution, from where
% the figures for the paper are generated, is already calculated and saved in the folder Data_solution.
% To plot the already solved and saved solution, run plot_all.m
% To newly calculate the solution or play with the code, run this file. Here more information:
% Parameters for the model are saved in the file params_stefan.m, 
% the macroscopic ODE is solved using ode15s with stefanG.m and fun_stefanG.m, 
% the microscopic ODE is solved using ode15s with stefanL.m and fun_stefanL.m (G for globa, L for local),
% the file event_melted_s.m is needed for the Maltab feature 'Events';
% Because it takes a while to calculate the solution, the function plot_help shows the progress.


clear all
close all
clc

global time;
global Tlocal;
global melted;
global pfw;
global pvw;
global dtS;


p=params_stefan;
time = 0;
melted = zeros(p.NG,1);
pfw = p.pfw0*ones(1,p.NG);
pvw = p.pvw0*ones(1,p.NG);
dtS = zeros(p.NG,1);
%Tlocal = cell(1,p.NG);
for k=1:p.NG
    Tlocal{k} = p.Tc*ones(1,p.NL);
end

solTemp = zeros(0,p.NG);
solR = zeros(0,p.NG);
solU = zeros(0,p.NG);
solSig = zeros(0,p.NG);
solSiw = zeros(0,p.NG);
solT = zeros(0,1);

%% ODE
timespan = [0 p.TimeMax];
Temp0 = p.TempInit*ones(p.NG,1);
r0 = p.r0*ones(p.NG,1);
U0 = p.U0*ones(p.NG,1);
sig0 = p.sig0*ones(p.NG,1);
siw0 = p.siw0*ones(p.NG,1);


y0 = [Temp0;r0;U0;sig0;siw0];

options = odeset('Events',@event_melted_s,'AbsTol',[p.eps*1e-5*ones(1,p.NG) p.eps*1e-5*ones(1,p.NG) p.eps^2*1e-10*ones(1,p.NG) p.eps*1e-5*ones(1,p.NG) p.eps*1e-5*ones(1,p.NG)]);


IE=0;
zahl=1;
TEvec = zeros(1,p.NG);
tic
while any(melted == 0)
    %if any of the ice bars is melted, an event happens, and the
    %calculation stops. Then the while-loop starts again.
    
    [T,Y,TE,YE,IE] = ode15s(@fun_stefanG,timespan,y0,options,p);
    
    
    indexT = find(T==TE);
    
    solTemp = [solTemp; Y(1:indexT,1:p.NG)];
    solR = [solR; Y(1:indexT,p.NG+1:2*p.NG)];
    solU = [solU; Y(1:indexT,2*p.NG+1:3*p.NG)];
    solSig = [solSig; Y(1:indexT,3*p.NG+1:4*p.NG)];
    solSiw = [solSiw;Y(1:indexT,4*p.NG+1:5*p.NG)];
    solT = [solT; T(1:indexT)];
    melted(IE) = 1;
    
    timespan = [T(indexT) p.TimeMax];
    
    y0 = YE;
    for k=1:p.NG
        if melted(k) == 1;
            y0(4*p.NG + k) = 0;
        end
    end
    disp('Geschmolzen:')
    IE
    TEvec(zahl) = length(solT);
    zahl=zahl+1;
end
toc


options = odeset('RelTol',1e-4);
[T,Y] = ode15s(@fun_stefanG,timespan,y0,options,p);

solTemp = [solTemp; Y(:,1:p.NG)];
solR = [solR; Y(:,p.NG+1:2*p.NG)];
solU = [solU; Y(:,2*p.NG+1:3*p.NG)];
solSig = [solSig; Y(:,3*p.NG+1:4*p.NG)];
solSiw = [solSiw;Y(:,4*p.NG+1:5*p.NG)];
solT_end = [solT; T];
pfw(2,:)=[];
pvw(2,:)=[];

j=8;
for i=1:j-1
    TEhelp(i) = length(solT) + floor(length(T)/j)*i;
end
TEvec = [TEvec TEhelp];
%% save result

save('solTemp.mat','solTemp') 
save('solR.mat','solR')
save('solU.mat','solU')
save('solSig.mat','solSig')
save('solSiw.mat','solSiw')
save('solT.mat','solT')
save('solT_end.mat','solT_end')
save('pfw.mat','pfw')
save('pvw.mat','pvw')
save('TEvec.mat','TEvec')



plot_thaw3(solT,solTemp,solR,solU,solSig,solSiw,pfw,pvw,p)
plot_thaw2(solT,solT_end,solTemp,solR,solU,solSig,solSiw,pfw,pvw)



