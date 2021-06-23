function hwtr=h_watergap(Tin,Tout,m,wgap,wc,Lc)
% the heat transfer coefficienct of laminar flow between two flat plates.
% reference: Anthony, L.H. and R.N. Maddox, Mass transfer: fundamentals andapplications. 1985: Prentice-Hall Englewood-Cliffs, NJ.
% Tin=350;
% Tout=370;
% m=0.01;
% wgap=0.01*1.6*0.95*2;
% wc=0.95;
% Lc=1.6;

Tavg=(Tin+Tout)/2;%average temperature
vw=exp(-3.7188+578.919/(-137.546+0.5*Tin+0.5*Tout))/1000;%viscosity of water
kw=(9.28516*10^(-7)*Tavg^3 - 1.06167*10^(-2)*Tavg^2 + 7.76041*Tavg - 7.87144*10^2)/1000;% water thermal conductivity
Prw=vw*4180/kw; % water Pr
Rew=(m/wgap/wc)*(2*wgap)/vw;% water Reynolds number in pipe; the hydraulic diameter is 4A/C=2wgap

%the length of the thermal entrance
Lentr=0.05*Rew*Prw*(2*wgap);
h_entrance=(7.54+0.03*(2*wgap/Lentr)*Rew*Prw/(1+0.016*((2*wgap/Lentr)*Rew*Prw)^(2/3)))*kw/(2*wgap);%heat transfer coefficnect in the entrance part
h_fulldev=7.54*kw/(2*wgap);%heat transfer coefficnect in the full-developed entrance part

if Lentr<Lc
    hwtr=h_entrance*(Lentr/Lc)+h_fulldev*((Lc-Lentr)/Lc);
else
    hwtr=h_entrance;
end
    

end