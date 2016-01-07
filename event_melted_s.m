function [value,isterminal,direction]=event_melted_s(~,y,p)
global melted

Sig = y(3*p.NG+1:4*p.NG);
Siw = y(4*p.NG+1:5*p.NG);

zahl = 1;
value = ones(p.NG,1);
for k=1:p.NG
    if melted(k) == 0
        value(k) = Siw(zahl)-Sig(zahl);
        zahl = zahl + 1;
    elseif melted(k) == 1
        value(k) = 1;
    end
end

isterminal=ones(p.NG,1);
direction=zeros(p.NG,1);