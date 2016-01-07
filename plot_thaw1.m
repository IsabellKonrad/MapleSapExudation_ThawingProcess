function plot_thaw1(solT,solTemp,solR,solU,solSig,solSiw,pfw,pvw)


p=params_stefan;
solTemp = [solTemp ones(size(solTemp,1),1)*p.TempOut];
langTemp = [p.langG p.LG]';
short_time = 1:4:length(solT);
short_time(end+1) = length(solT);



%% Plots
% Please uncomment requested issue


%% Temperature
figure
for j=short_time
    plot(langTemp,solTemp(j,:),'k','LineWidth',3);
    axis([0 p.LG p.Tc p.TempOut])
    title(sprintf('Temperature,   time %i',solT(j)),'FontSize',16);
    xlabel('Radius of tree trunk in m','FontSize',14)
    ylabel('T in K','FontSize',14)
    pause(0.03)
end


%% Sig and Siw
for k=1:p.NG
    indexhelp = find(solSiw(:,k) == 0,1,'first');
    solSiw(indexhelp:end,k) = solSig(indexhelp:end,k);
end
figure
maxSiw = max(max(solSiw));
minSig = min(min(solSig));
for j=short_time
    plot(p.langG,solSig(j,:),'b','LineWidth',3);
    hold on
    plot(p.langG,solSiw(j,:),'k','LineWidth',3);
    axis([0 p.LG minSig maxSiw])
    title(sprintf('S_{gi} and S_{iw},   time %i',solT(j)),'FontSize',16);
    xlabel('Radius of tree trunk in m','FontSize',14)
    ylabel('Radius in m','FontSize',14)
    hold off
    pause(0.02)
end



%% Radius R
figure
maxR = max(max(solR));
minR = min(min(solR));
for j=short_time
    plot(p.langG,solR(j,:),'k','LineWidth',3);
    axis([0 p.LG minR maxR])
    title(sprintf('Radius R,   time %i',solT(j)),'FontSize',14);
    xlabel('Radius of tree trunk in m','FontSize',14)
    ylabel('Radius in m','FontSize',14)
    pause(0.02)
end


%% Water volume U
figure
maxU = max(max(solU));
minU = min(min(solU));
for j=short_time
    plot(p.langG,solU(j,:),'k','LineWidth',3);
    axis([0 p.LG minU maxU])
    title(sprintf('Water Volume U,   time %i',solT(j)),'FontSize',16);
    xlabel('Radius of tree trunk in m','FontSize',14)
    ylabel('Volume in m^3','FontSize',14)
    pause(0.02)
end


%% Water Pressure in Fiber
figure
maxPfw = max(max(pfw));
minPfw = min(min(pfw(2:end,:)));
for j=short_time
    plot(p.langG,pfw(j,:),'k','LineWidth',3);
    axis([0 p.LG minPfw maxPfw])
    title(sprintf('Presure water in Fiber,   time %f',solT(j)),'FontSize',13,'FontWeight','bold');
    pause(0.02)
end


% Water pressure in Vessel
figure
maxPvw = max(max(pvw));
minPvw = min(min(pvw(2:end,:)));
for j=short_time
    plot(p.langG,pvw(j,:),'k','LineWidth',3);
    axis([0 p.LG minPvw maxPvw])
    title(sprintf('Presure water in Vessel,   time %i',solT(j)),'FontSize',16);
    xlabel('Radius of tree trunk in m','FontSize',14)
    ylabel('Pressure p_w^v','FontSize',14)
    pause(0.02)
end

% movie2avi(film1, 'film_temp.avi', 'compression', 'None');
% movie2avi(film2, 'film_sigsiw.avi', 'compression', 'None');
% movie2avi(film3, 'film_r.avi', 'compression', 'None');
% movie2avi(film4, 'film_U.avi', 'compression', 'None');
% movie2avi(film5, 'film_pfw.avi', 'compression', 'None');
% movie2avi(film6, 'film_pvw.avi', 'compression', 'None');

end