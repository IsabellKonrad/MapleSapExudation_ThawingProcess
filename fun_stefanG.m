function dy = fun_stefanG(t,y,p)
global time;
global Tlocal;
global melted;
global pfw;
global pvw;
global dtS;


Temp=y(1:p.NG);
R = y(p.NG+1:2*p.NG);
U = y(2*p.NG+1:3*p.NG);
Sig = y(3*p.NG+1:4*p.NG);
Siw=y(4*p.NG+1:end);


%% Delete the values, if ode15s reverses time.
wrong = find(t<=time);
if( ~isempty(wrong) && t>0)
    time(wrong)=[];
    pfw(wrong,:)=[];
    pvw(wrong,:)=[];
    try
        for k=1:p.NG
            Tlocal{k}(wrong,:) = [];
        end
    catch
        warning('MATLAB:FUN_Stefan','Matlab tries to delete non-existing Tlocal values')
    end
end

%% Parameters for the next step
pressure_f_w = p.pfg0*p.sig0^2./Sig.^2 - 2*p.sigmaw./Sig;
Vvg = pi*p.Lv*R.^2;
rhovg = p.rhovg0*p.Vvg0./(Vvg + p.H*(p.Vv - Vvg));
pvg = rhovg.*Temp*p.R/p.Mg;
pressure_v_w = pvg - 2*p.sigmas./R;
dU = -p.K*p.A/(p.Nf*p.rhow*p.g*p.W)*(pressure_v_w - pressure_f_w - p.R*Temp*p.cs);

%% Calculate the local temperature
[dxTout,dxTin] = stefanL(t,Temp,Siw,p);


%% Calculation for the geometry in the reference cell
Gamma1 = p.Dw*2*pi*p.r2*ones(p.NG,1)/p.eps;
Y1 = p.eps-Siw.^2*pi/p.eps;
%consider heat diffusion in gas:
PT = p.eps*p.function_difftensor(Siw,R);
%PT = p.eps*p.function_PT(Siw);

index = find(melted == 1);
Y1(index) = 1;
Gamma1(index) = 0;
%consider heat diffusion in gas:
PT(index) = p.function_difftensor_melted(Sig(index),R(index));
%PT(index) = p.function_PT(Siw(index)*0);
 
 %% Derivative of temperature and ice bar
 dSiw = zeros(p.NG,1);
 dSig = zeros(p.NG,1);
 dR = zeros(p.NG,1);
 
 dTemp = ((p.AG*Temp + p.bG).*PT - Gamma1.*dxTout)./Y1/p.cw;
 
 for k = 1:p.NG
     if melted(k) == 0
         dSiw(k) =  -p.kw/(p.lambda*p.rhow)*dxTin(k) + dU(k)/(2*pi*Siw(k)*p.Lf);
         dSig(k) = dU(k)*p.rhow/(2*p.rhoi*Sig(k)*pi*p.Lf) + 2*(p.rhoi - p.rhow)*Siw(k)*dSiw(k)/(2*p.rhoi*Sig(k));
         dR(k) = -dU(k)*p.Nf/(2*pi*R(k)*p.Lf);
         if (Siw(k) >= p.Rf && dSiw(k)>0)
             dSiw(k) = 0;
             dU(k) = p.kw/(p.lambda*p.rhow)*dxTin(k)*(2*pi*Siw(k)*p.Lf);
             dSig(k) = dU(k)*p.rhow/(2*p.rhoi*Sig(k)*pi*p.Lf);
             dR(k) = -dU(k)*p.Nf/(2*pi*R(k)*p.Lf);
         end
     elseif melted(k) == 1
         dSiw(k) = 0;
         dSig(k) = dU(k)/(2*pi*Sig(k)*p.Lf);
         dR(k) = -dU(k)*p.Nf/(2*pi*R(k)*p.Lf);
     end
 end
 
dtS = dSiw;

dy = [dTemp;dR;dU;dSig;dSiw];

%% Remember next time step
if (t>0)
    time(end+1) = t;
    pfw(end+1,:) = pressure_f_w;
    pvw(end+1,:) = pressure_v_w;
end
end
