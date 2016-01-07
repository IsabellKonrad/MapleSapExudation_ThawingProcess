function plot_help(solT,time_vec,solTemp,solSig,solSiw)

p=params_stefan;
solTemp = [solTemp ones(size(solTemp,1),1)*p.TempOut];
langTemp = [p.langG p.LG]';
indexShort_time = time_vec;
indexShort_time(end+1) = length(solT);


plot(langTemp,solTemp(end,:),'k','LineWidth',3);
    axis([0 p.LG p.Tc p.TempOut])
    title(sprintf('Temperature,   time %i',solT(end)),'FontSize',16);
    xlabel('Radius of tree trunk in m','FontSize',14)
    ylabel('T in K','FontSize',14)


for k=1:p.NG
    indexhelp = find(solSiw(:,k) == 0,1,'first');
    solSiw(indexhelp:end,k) = solSig(indexhelp:end,k);
end

figure
maxSiw = max(max(solSiw));
minSig = min(min(solSig));
axis([0 p.LG minSig maxSiw])
title(sprintf('S_{gi} and S_{iw}'),'FontSize',16);
xlabel('Radius of tree trunk in m','FontSize',14)
ylabel('Radius in m','FontSize',14)
hold on
for j=indexShort_time
    plot(p.langG,solSig(j,:),'b','LineWidth',3);
    plot(p.langG,solSiw(j,:),'k','LineWidth',3);
end
end