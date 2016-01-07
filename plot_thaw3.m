function plot_thaw3(solT,solTemp,solR,solU,solSig,solSiw,pfw,pvw)

p=params_stefan;


%% Plotting
set( 0, 'DefaultTextFontName', 'times' );
set( 0, 'DefaultAxesFontName', 'times' );
set(0, 'DefaultFigurePosition', [150 150 900 500] );

set( 0, 'DefaultTextFontSize',  24 );
set( 0, 'DefaultAxesFontSize',  24 );
set( 0, 'DefaultLineLineWidth', 2 );   % default is 1
set( 0, 'DefaultLineMarkerSize', 8 );   % default is 6



solTemp = [solTemp ones(size(solTemp,1),1)*p.TempOut];
langTemp = [p.langG p.LG]';
%zahl=1;
% for k = 1:length(TEvec)
%     if 0==mod(k-1,3)
%         TEvec_kurz(zahl) = TEvec(k);
%         zahl=zahl+1;
%     end
% end
% TEvec_kurz(end) = TEvec(end);
% TEvec = 0;
% TEvec = TEvec_kurz;
% TEvec=[TEvec 851 856 859];

yy = size(solT,1);
vec = [2:yy/5:yy yy];
zz = floor(vec);
stylevec = cell(1,6);
stylevec{1} = '-.';
stylevec{2} = '-';
stylevec{3} = '--';
stylevec{4} = '-.';
stylevec{5} = '-';
stylevec{6} = '--';
%grayvec = [[0 0 0]; [0.12 0.12 0.12]; [0.24 0.24 0.24]; [0.36 0.36 0.36]; [0.48 0.48 0.48] ; [0.6 0.6 0.6]];
colorvec = [[220 0 0]; [0 130 0]; [0 0 150]; [170 0 170]; [200 100 0] ; [0 100 200]]/255;


%% Temperature
fig1 = figure;
hold on
axis([0 p.LG p.Tc p.TempOut])
zahl = 1;
for j=zz
    plot(langTemp,solTemp(j,:),stylevec{zahl},'Color',colorvec(zahl,:));
    zahl = zahl + 1;
end
xlabel('$x\ [m]$','Interpreter','latex');
ylabel('$T_1\ [{}^\circ K]$','Interpreter','latex');
leg1 = legend(sprintf('%0.1f h',solT(zz(1))/3600),...
    sprintf('%0.1f h',solT(zz(2))/3600),...
    sprintf('%0.1f h',solT(zz(3))/3600),...
    sprintf('%0.1f h',solT(zz(4))/3600),...
    sprintf('%0.1f h',solT(zz(5))/3600),...
    sprintf('%0.1f h',solT(zz(6))/3600),'Location','NorthWest');
set(leg1, 'Box', 'off')

line([0.23 0.11],[278 278],'Color','k');
line([0.11  0.115],[278 278.5],'Color','k'); 
line([0.11  0.115],[278 277.5],'Color','k'); 
text('String','increasing time','Position',[0.13 278.5])

print(fig1,'-depsc','Temp_hom_thaw.eps')



%% Sig and Siw
for k=1:p.NG
    indexhelp = find(solSiw(:,k) == 0,1,'first');
    solSiw(indexhelp:end,k) = solSig(indexhelp:end,k);
end
fig2 = figure;
maxSiw = max(max(solSiw(1:zz(end))))*1.01;
minSig = min(min(solSig(1:zz(end))))*0.99;
axis([0 p.LG minSig*1e6 maxSiw*1e6])
hold on
zahl=1;
for j=zz
    plot(p.langG,solSig(j,:)*1e6,stylevec{zahl},'Color',colorvec(zahl,:));
    zahl = zahl+1;
end
xlabel('$x\ [m]$','Interpreter','latex');
ylabel('$s\ [\mu m]$','Interpreter','latex');

zahl=1;
for j=zz
    plot(p.langG,solSiw(j,:)*1e6,stylevec{zahl},'Color',colorvec(zahl,:));
    zahl = zahl+1;
end
text('Interpreter','late','String','$s_{iw}$','Position',[0.21 3.45])
text('Interpreter','late','String','$s_{gi}$','Position',[0.21 2.53])
hold off
print(fig2,'-depsc','s_hom_thaw.eps')


%% Radius R
fig3 = figure;
maxR = max(max(solR(1:zz(end),:)))*1.01;
minR = min(min(solR(1:zz(end),:)))*0.99;
axis([0 p.LG minR*1e6 maxR*1e6])
hold on
zahl = 1;
for j=zz
    plot(p.langG,solR(j,:)*1e6,stylevec{zahl},'Color',colorvec(zahl,:));
    zahl = zahl+1;
end
xlabel('$x\ [m]$','Interpreter','latex');
ylabel('$r\ [\mu m]$','Interpreter','latex');
print(fig3,'-depsc','r_hom_thaw.eps')



%% Water volume U
fig4 = figure;
maxU = max(max(solU(1:zz(end),:)))*1.01;
minU = min(min(solU(1:zz(end),:)))*0.99;
axis([0 p.LG minU*1e18 maxU*1e18])
hold on
zahl = 1;
for j=zz
    plot(p.langG,solU(j,:)*1e18,stylevec{zahl},'Color',colorvec(zahl,:));
    zahl = zahl+1;
end
xlabel('$x\ [m]$','Interpreter','latex');
ylabel('$U\ [\mu m^3]$','Interpreter','latex');
print(fig4,'-depsc','U_hom_thaw.eps')



%% Water Pressure in Fiber
fig5 = figure;
maxPfw = max(max(pfw(1:zz(end),:)))*1.01;
minPfw = min(min(pfw(2:zz(end),:)))*0.99;
axis([0 p.LG minPfw/1e3 maxPfw/1e3])
hold on
zahl = 1;
for j=zz
    plot(p.langG,pfw(j,:)/1e3,stylevec{zahl},'Color',colorvec(zahl,:));
    zahl=zahl+1;
end
xlabel('$x\ [m]$','Interpreter','latex');
ylabel('$p_w^f\ [kPa]$','Interpreter','latex');
print(fig5,'-depsc','pfw_hom_thaw.eps')



%% Water pressure in Vessel
fig6 = figure;
maxPvw = max(max(pvw(1:zz(end),:)))*1.01;
minPvw = min(min(pvw(2:zz(end),:)))*0.99;
axis([0 p.LG minPvw/1e3 maxPvw/1e3])
zahl = 1;
hold on
for j=zz
    plot(p.langG,pvw(j,:)/1e3,stylevec{zahl},'Color',colorvec(zahl,:));
    zahl = zahl+1;
end
xlabel('$x\ [m]$','Interpreter','latex');
ylabel('$p_w^v\ [kPa]$','Interpreter','latex');
print(fig6,'-depsc','pvw_hom_thaw.eps')

end