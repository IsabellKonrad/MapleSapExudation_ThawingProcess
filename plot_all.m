% Run this file to plot the already solved solution saved in Data_solution. 
% The solution in Data_solution is used to generate the figures for the Paper
% Homogenization of the Stefan problem, with application to maple sap exudation.
% This file calls the three plotting-helper-functions plot_thaw1, plot_thaw2 and plot_thaw 3
% plot_thaw1 makes short videos from the temperature, the radius of the ice in the fiber Siw and
% the radius of the gas in the fiber Sig, the radius of the gas in the vessel R, the water volume
% moving from the fiber into the vessel U and the water pressure in fiber pfw and vessel pvw.
% plot_thaw2 makes figures with x-axis the time and y-axis the temperature, Siw and Sig, 
% R, U and pfw and pvw.
% plot_thaw3 makes figures with x-axis the radius of the tree and y-axis the temperature, Siw and Sig,
% R, U, pfw and pvw.



load('Data_solution/pfw.mat')
load('Data_solution/pvw.mat')
load('Data_solution/solR.mat')
load('Data_solution/solSig.mat')
load('Data_solution/solSiw.mat')
load('Data_solution/solT.mat')
load('Data_solution/solT_end.mat')
load('Data_solution/solTemp.mat')
load('Data_solution/solU.mat')
load('Data_solution/TEvec.mat')

plot_thaw1(solT,solTemp,solR,solU,solSig,solSiw,pfw,pvw)
plot_thaw2(solT,solT_end,solTemp,solR,solU,solSig,solSiw,pfw,pvw)
plot_thaw3(solT,solTemp,solR,solU,solSig,solSiw,pfw,pvw)
