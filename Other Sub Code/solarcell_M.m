function CellCal=solarcell(Ac,Gpv,G,cell)
%Tpv: cell temperature, K
%Ac: PV total area
%Gpv:Spectrum arrived the PV
%G: actual solar irrandiance of the ambient
%cell:cell type



% test data for this function
% clc;clear;
% Tpv=273.15+25;
% Ac=0.95*1.6;
% Gpv=xlsread('data', 'ambient conditions', 'C2:C2003'); %the distribution of AM1.5
% G=1000;
% CellCal=solarcell(Ac,Gpv,G,'silicon')
% CellCal=solarcell(Ac,Gpv,G,'CdTe')


%references:
%(1)Fan, J. C. (1986). Theoretical temperature dependence of solar cell parameters. Solar cells, 17(2-3), 309-315.
%(2)Otanicar, T., et al., Parametric analysis of a coupled photovoltaic/thermal concentrating solar collector for electricity generation. Journal of Applied Physics, 2010. 108(11): p. 114907.

Tpv0=30.79+273.15;
kB=1.38064852E-23; %Boltzman constant
e=1.602176634E-19; % charge of electron,constant

if G~=0
    
    if strcmp(cell,'silicon')==1
        k1=0.03;%0.63 %empirical parameter, Fan, J. C. (1986). Theoretical temperature dependence of solar cell parameters. Solar cells, 17(2-3), 309-315.
        b=1.25; %empirical parameter
        n=2.98; %0.71%empirical parameter 

    Ebg=1.11*e; % bandgap energy, eV
    A1=1; % cell ideality factor,n=1 for the single p-n junction, A1 should be 1~2
    Tcoeff=0.0045;%0.0045
    
    SR=xlsread('data', 'solar cells', 'I3:I2004'); % Si, SunPowe C60 Solar Cell,cell efficiency 25.1%
    %simulation result:25.3%
    end
    
    
    if strcmp(cell,'CdTe')==1  
        k1=0.03; %empirical parameter
        b=1.15; %empirical parameter
        n=0.98; %empirical parameter 

    Ebg=1.44*e; % bandgap energy, eV
    A1=1; % cell ideality factor,n=1 for the single p-n junction, A1 should be 1~2
    Tcoeff=0.0025;
    
    SR=xlsread('data', 'solar cells', 'D3:D2004'); % CdTe, First Solar 4 sereries,cell efficiency 22.1%
    %simulation result:22.2%
    end
    
 
    
    N=size(SR);
    J0=k1*Tpv0^(3/n)*exp(-Ebg/(b*kB*Tpv0));
    Jsc=0;
    lamadaG=xlsread('data', 'ambient conditions', 'B2:B2004');

        for i=1:N(1)
           Jsc=Jsc+Ac*Gpv(i)*SR(i)*(lamadaG(i+1)-lamadaG(i)); 
        end
    
    Voc=A1*kB*Tpv0/e*log(Jsc/J0+1);
    voc=Voc/(A1*kB*Tpv0/e);
    FF=(voc-log(voc+0.72))/(voc+1);
    Effm=Voc*Jsc*FF/(G*Ac);%effm is the module efficiency at 25 oC
    
     Gpv1=0;%the energy of the solar arrived pv
    
        for j=1:N(1)
            Gpv1=Gpv1+Gpv(j)*(lamadaG(j+1)-lamadaG(j));
        end

    
    %Gabs=Gpv1-Voc*Jsc*FF*(1-Tcoeff*(Tpv-273.15-25))/Ac;%the waste heat in pv, W/m2
  

    CellCal(1)=Effm;
    CellCal(2)=FF;
    CellCal(3)=Voc;
    CellCal(4)=Jsc;
    CellCal(5)=Gpv1;
    CellCal(6)=Tcoeff;

else 
    CellCal(1)=0;
    CellCal(2)=0;
    CellCal(3)=0;
    CellCal(4)=0;
    CellCal(5)=0;
    CellCal(6)=0;
end

end