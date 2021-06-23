function [K,v,B,Pr]=AirProperties(Tprop)
%air properties as function of the temperatures
% %air properties REF Parametric Studies of Top Heat Loss Coefficiency of
% %Double glazed flat plate solar collector, Bisen,Dass,Jain 2011
 K=-3e-8*Tprop^2 + 1e-4*Tprop -4e-5;
 v=(9e-5*Tprop^2 + 0.040*Tprop -4.17)*1e-6;
 B=1/Tprop;
 Pr=1.057-0.06*log(Tprop);
end
