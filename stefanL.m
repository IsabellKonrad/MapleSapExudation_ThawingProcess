function [dxTout,dxTin] = stefanL(t,Temp,Siw,p)
global time;
global Tlocal;
global melted;
global dtS;

dxTout = zeros(p.NG,1);
dxTin = zeros(length(Siw),1);


for k = 1:p.NG
    if (time(end)<t && melted(k) == 0)
        % moving mesh method used
        
        b = p.r2;
        vector = dtS(k)*(p.langL-1)/(b-Siw(k)) + p.Dw./(p.langL*(b-Siw(k))+Siw(k))/(b-Siw(k));
        A4 = diag(vector(2:p.NL),1) - diag(vector(3:p.NL+1),-1);
        AL = p.Dw/(b-Siw(k))^2/p.hL^2*p.A3 + A4/2/p.hL;
        bL = [p.Tc*(p.Dw/(b-Siw(k))^2/p.hL^2 - vector(2)/2/p.hL);zeros(p.NL-2,1);Temp(k)*(p.Dw/(b-Siw(k))^2/p.hL^2 + vector(end-1)/2/p.hL)];
   
        %ODE
        start_u = Tlocal{k}(end,:);
        timespan = [time(end) t];
        options = odeset('AbsTol',p.eps*1e-3*ones(1,p.NL));
        [~,u] = ode15s(@fun_stefanL,timespan,start_u,options,AL,bL,p);
        
        Tlocal{k}(end+1,:) = u(end,:);
        dxTout(k) = (Temp(k) - u(end,end))/p.hL/(b-Siw(k));
        dxTin(k) = (u(end,1) - p.Tc)/p.hL/(b-Siw(k));
 
    elseif (t>0)
        Tlocal{k}(end+1,:) = Tlocal{k}(end,:);
    end
end

end