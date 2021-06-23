function h_gap = h_enclosure(Ttop,Tbot,gap,beta)
%calculation of heat transfer coefficienct in enclosed surface with natural
%convection 
Tmean=(Ttop+Tbot)/2;
[K,v,B,Pr]=AirProperties(Tmean);
%[K,v,B,Pr]=ArgonProperties(Tmean);

 if (Ttop==Tbot)
     h_gap=0+1e-9;     
 else
 Ra=9.8*abs(Ttop-Tbot)*gap^3*Pr/v^2/Tmean;
 
 A1=1-1708/Ra/cos(beta);
 B1=(Ra*cos(beta)/5830)^(1/3)-1;
 C1=1-1708*(sin(1.8*beta))^1.6/Ra/cos(beta);
 if A1<0
    A=0;
 else
    A=A1;
 end
 
if B1<0
      Bb=0;
else
      Bb=B1;
end

if C1<0
      C=0;
else
      C=C1;
end
Nu=1+1.44*A*C+Bb; %reference "Natural convection in enclosed spaces a review of application to solar energy collection"

Tmeangap=(Ttop+Tbot)/2;
[KNu,vNu,BNu,PrNu]=AirProperties(Tmeangap);
h_gap=Nu*KNu/gap;
end
end

