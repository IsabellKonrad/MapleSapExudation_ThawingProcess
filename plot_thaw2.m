function plot_thaw2(solT,solT_end,solTemp,solR,solU,solSig,solSiw,pfw,pvw)

p=params_stefan;
%allDiscP = 1:p.NG;
%FourDiscP = [10 20 30 40 50];
FourDiscP = 30;
bis = size(solT,1);

%% Plotting
set( 0, 'DefaultTextFontName', 'times' );
set( 0, 'DefaultAxesFontName', 'times' );
set(0, 'DefaultFigurePosition', [150 150 800 500] );

set( 0, 'DefaultTextFontSize',  24 );
set( 0, 'DefaultAxesFontSize',  24 );
set( 0, 'DefaultLineLineWidth', 2 );   % default is 1
set( 0, 'DefaultLineMarkerSize', 8 );   % default is 6
set( 0, 'DefaultLineColor', 'k');


%% Temperature Temp
fig1 = figure;
maxTemp = max(solTemp(1:bis,FourDiscP))*1.001;
minTemp = min(solTemp(1:bis,FourDiscP))*0.999;
for j=FourDiscP
    plot(solT_end(1:bis)/3600,solTemp(1:bis,j),'Color',[0 130 0]/255);
end
axis([0 solT_end(bis)/3600 minTemp maxTemp]);
xlabel('$t\ [h]$','Interpreter','latex')
ylabel('$T_1\ [{}^\circ K]$','Interpreter','latex');
print(fig1,'-depsc','Temp_30.eps')



%% Siw and Sig
for k=1:p.NG
    indexhelp = find(solSiw(:,k) == 0,1,'first');
    solSiw(indexhelp:end,k) = solSig(indexhelp:end,k);
end

fig2 =figure;
maxS = max(solSiw(1:bis,FourDiscP))*1.03;
minS = min(solSig(1:bis,FourDiscP))*0.97;
for j=FourDiscP
    plot(solT_end(1:bis)/3600,solSig(1:bis,j)*1e6,'-','Color',[0 130 0]/255);
    hold on
    plot(solT_end(1:bis)/3600,solSiw(1:bis,j)*1e6,'-.','Color',[0 100 200]/255);
end
axis([0 solT_end(bis)/3600 minS*1e6 maxS*1e6]);
leg1 = legend('$s_{gi}$', '$s_{iw}$','Location','NorthEast');
set(leg1,'Interpreter','latex')
set(leg1, 'Box', 'off')
xlabel('$t\ [h]$','Interpreter','latex')
ylabel('$[\mu m]$','Interpreter','latex')
hold off
print(fig2,'-depsc','s_30.eps')



% Radius R
fig3 = figure;
for j=FourDiscP
    plot(solT_end(1:bis)/3600,solR(1:bis,j)*1e6,'k');
    xlabel('$t\ [h]$','Interpreter','latex')
    ylabel('$r\ [\mu m]$','Interpreter','latex')
end
print(fig3,'-depsc','r_30.eps')
 
 
 
 % Water volume U
fig4 = figure;
for j=FourDiscP
    plot(solT_end(1:bis)/3600,solU(1:bis,j)*1e18,'k');
    xlabel('$t\ [h]$','Interpreter','latex')
    ylabel('$U\ [\mu m^3]$','Interpreter','latex')
end
print(fig4,'-depsc','U_30.eps')




%% Water pressure in fiber  
fig5 = figure;
maxP = max(pvw(1:bis,FourDiscP))*1.05;
minP = min(pfw(1:bis,FourDiscP))*0.95;
for j=FourDiscP
    plot(solT_end(1:bis)/3600,pfw(1:bis,j)/1e3,'-','Color',[0 130 0]/255);
    hold on
    plot(solT_end(1:bis)/3600,pvw(1:bis,j)/1e3,'-.','Color',[0 100 200]/255);
end
axis([0 solT_end(bis)/3600 minP/1e3 maxP/1e3]);
leg1 = legend('$p_w^f$', '$p_w^v$','Location','East');
set(leg1,'Interpreter','latex')
set(leg1, 'Box', 'off')
xlabel('$t\ [h]$','Interpreter','latex')
ylabel('$[kPa]$','Interpreter','latex')
print(fig5,'-depsc','p_30.eps')
end
