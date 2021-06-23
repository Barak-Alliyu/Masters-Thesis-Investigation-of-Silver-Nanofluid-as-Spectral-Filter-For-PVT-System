 %Spectral splitting PVT steady-state code: a filter covered on a comercial PV (3/7/2019)
%The structure of the S-PVT
   
%=========CHANNEL TOP GLASS==========   T1
%           HEAT TRANSFER FLUID         T2:inlet  T3:outlet
%------------FILTER TOP GLASS--------   T4
%            FILTER FLUID               T5
%------------FILTER BOTTOM GLASS-----   T6
%          HEAT TRANSFER FLUID          T2:inlet  T9:outlet
%========CHANNEL BOTTOM GLASS========   T7
%             AIR GAP                    
%----------------PV------------------   T8

%clc;clear;
%Geometry parameters
Lc=0.116; % length - Aperture length, m                           
Wc=0.116; % Aperture width, m
Nc=1; % number of collectors
Ac=Lc*Wc*Nc; %Lc*Wc; % aperture area, m2
Asc=0.003375; % area of solar cell
Wgap1=0.01;%The gap thickness between channel top glass and the filter top glass
Wgap2=0.01;%The gap thickness between filter top glass and the filter bottom glass
Wgap3=0.01;%The gap thickness between the filter bottom glass and channel bottom glass
Wgap4=0.005;%The gap thickness between the channel bottom glass and the PV

%Tfl1=filter top glass
%Tfl2=filter bottom glass

cell='silicon';
FLdFlt='valvoline';% insert appropriate fluid
kfl=0.1384;%thermal conductivity of water
GG=[200 400 600 800 1000];
Ta=308.0;%300;%Ambient temperature K
Twtin=298.0;%Heat transfer fluid inlet temperature
mwt=[0.0079 0.008 0.009 0.01 0.011 0.012 0.013 0.014 0.015]*Ac;
%0.003715814507 0.007431629013 0.0111474]*Ac;%mass flow rate of the channel fd
% 0.0005 0.001 0.002 0.003 0.004 0.005 0.006



Wflg=0.0007;%thickness of the filter top and bottom glass
Wgc=0.0038;%thickness of the channel top and bottom glass


%Optical properties
%select working fluid filter
%layer:pure_water_10mm,anti_reflective_thin_glass_3.8mm, anti_reflective_thin_glass_0.7mm, nanofluid_1, nanofluid_2, nanofluid_3, nanofluid_4, nanofluid_5, valvoline, glycerol, glycerol_water_mixture, ideal_filter_glass_Si_300-1200


Emgc=0.9; % emissivity of the channel glass
[Trgc,Relgc,Abgc]=select_filter_M('anti_reflective_thin_glass_3.8mm');

   
%Channel fluid
[Trwt,Relwt,Abwt]=select_filter_M('pure_water_10mm');

%Filter glass
Emgfl=0.9; % emissivity of the filter glass (top & bottom)
[Trflg,Relflg,Abflg]=select_filter_M('anti_reflective_thin_glass_0.7mm');

%working fluid in the filter
if strcmp(FLdFlt,'nanofluid_1')==1
    [Trflf,Relflf,Abflf]=select_filter_M('nanofluid_1');
end

if strcmp(FLdFlt,'nanofluid_2')==1
    [Trflf,Relflf,Abflf]=select_filter_M('nanofluid_2');
end

if strcmp(FLdFlt,'nanofluid_3')==1
    [Trflf,Relflf,Abflf]=select_filter_M('nanofluid_3');
end


if strcmp(FLdFlt,'valvoline')==1
    [Trflf,Relflf,Abflf]=select_filter_M('valvoline');
end

if strcmp(FLdFlt,'glycerol')==1
    [Trflf,Relflf,Abflf]=select_filter_M('glycerol');
end

if strcmp(FLdFlt,'pure_water_10mm')==1
    [Trflf,Relflf,Abflf]=select_filter_M('pure_water_10mm');
end

Abpv=0.88; % 0.93 absorptivity of pv
Empv=0.75*ones(2002,1); % emmisivity of pv


%Ambient conditions and physical properties
%G=1000;%Solar intensity W/m2

AM=xlsread('data_M', 'ambient conditions', 'C2:C2003'); %the distribution of AM1.5
GAM=1000.37;%The intensity of the standard AM1.5

Vwind=0;%1;%Wind speed m/s

Tsky=0.0552*Ta^1.5; % sky temp

%mwt=[0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.015 0.02]*Ac;%mass flow rate of the channel fluid
% mwt=0.006193024178*Ac;
%h_cp=385/0.01;%thermal conductivity of copper/heat transfer
%T_cp=298;
cpwt=4180;%heat capacity of water
Tss=329;%temperature at solar simulator
Emss=0.35;%Emmisivity of solar simulator

Twt=Twtin;
%Emcp=0.07;
%****
%====================CALCULATION=============================
Nm=size(mwt);
NG=size(GG);


for j=5%j=1:NG(2)
G=GG(j);    
%Spectrum spliting on each layer
G0=G/GAM*AM;%original incident full-spectrum solar

[G1pass,G1rel,G1abs]=Gfilter(G0,Trgc,Relgc,Abgc);%layer1:top channel glass,Spectrum after the top channel glass
[G2pass,G2rel,G2abs]=Gfilter(G1pass,Trwt,Relwt,Abwt);%layer2:channel fluid(water_10mm),Spectrum after the top channel fluid
[G3pass,G3rel,G3abs]=Gfilter(G2pass,Trflg,Relflg,Abflg);%layer3:top filter glass,Spectrum after the top filter glass
[G4pass,G4rel,G4abs]=Gfilter(G3pass,Trflf,Relflf,Abflf);%layer4:fluid filter,Spectrum after fluid filter
[G5pass,G5rel,G5abs]=Gfilter(G4pass,Trflg,Relflg,Abflg);%layer5:bottom filter glass,Spectrum after the bottom filter glass
[G6pass,G6rel,G6abs]=Gfilter(G5pass,Trwt,Relwt,Abwt);%layer 6: channel fluid (water_10mm), Spectrum after bottom channel fluid
[G7pass,G7rel,G7abs]=Gfilter(G6pass,Trgc,Relgc,Abgc);%layer 7:bottom channel glass, Spectrum after bottom channel glass: G7pass is the spectrum arrived the pv


%each layer will absorb the reflective spectrum
[G7pass1,G7rel1,G7abs1]=Gfilter(G7pass*(1-Abpv),Trgc,Relgc,Abgc);%reflection reaching the bottom channel glass()
G7abs=G7abs+G7abs1;

[G6pass1,G6rel1,G6abs1]=Gfilter(G7rel+G7pass1,Trwt,Relwt,Abwt);%reflection reaching the bottom channel water
G6abs=G6abs+G6abs1;

[G5pass1,G5rel1,G5abs1]=Gfilter(G6rel+G6pass1,Trflg,Relflg,Abflg);%reflection reaching the bottom filter glass
G5abs=G5abs+G5abs1;

[G4pass1,G4rel1,G4abs1]=Gfilter(G5rel+G5pass1,Trflf,Relflf,Abflf);%reflection reaching the filter fluid
G4abs=G4abs+G4abs1;

[G3pass1,G3rel1,G3abs1]=Gfilter(G4rel+G4pass1,Trflf,Relflf,Abflf);%reflection reaching the top filter glass
G3abs=G3abs+G3abs1;

[G2pass1,G2rel1,G2abs1]=Gfilter(G3rel+G3pass1,Trwt,Relwt,Abwt);%reflection reaching the top channel water
G2abs=G2abs+G2abs1;

[G1pass1,G1rel1,G1abs1]=Gfilter(G2rel+G2pass1,Trgc,Relgc,Abgc);%reflection reaching the top channel glass
G1abs=G1abs+G1abs1;


% the performance of PV cell
Solarcell=solarcell(Asc,G7pass*Abpv,G,cell);%()

for i=1:Nm(2)
    
%simulation initials and settings
Asol=zeros(9,9); % 11 Equations' left-hand coefficients
Bsol=zeros(9,1); % 11 equations' right-hand coefficients
T=zeros(1,9)+300;% T1...T9:Please see the structure of the S-PVT
X=zeros(9,1)+300;% solution of each iteration
Ttemp=9999;
N=1;
error=1e-4;

while (max(abs(X.'-Ttemp))>error && N<100) 
    Ttemp=T;
% heat transfer coefficients calculations
[Kair,vair,Bair,Prair]=AirProperties(0.5*T(7)+0.5*T(8)); %air gap between channel bottom glass and pv
hr_gc2ss=0.01*h_radiation(Emgc,Emss,T(1),Tss); % radiation from channel top glass to solar simulator %0.06211 View factor from literature for 20cm distance and 16cm by 10cm
hwind=2*4.5+2.9*Vwind; %8.3+2.2*vwind(i-1); % 4.5+2.9*vwind, convection heat transfer due to wind
hr_gc2pv=0.01*h_radiation(Empv(1),Emgc,T(8),T(7)); % radiation between channel bottom glass and pv
hc_gc2pv=1/(2/h_enclosure(T(7),T(8),Wgap4,0)+Wgap4/Kair); % convection and conduction between channel bottom glass and pv through air gap ()Wgap1
hr_gc22ss=0.01*h_radiation(Emgc,Emss,T(7),Tss); % ideal radiation from channel bottom glass to sky
%hc_gc2am=1/(2/h_enclosure(T(7),T_cp,Wgap4,0)+Wgap4/Kair); % convection and conduction between channel bottom glass and copper
hc_wt1=h_watergap(T(2),T(3),0.1*mwt(i),Wgap1,Wc,Lc); % heat transfer coefficient in top water channel 
hc_wt2=h_watergap(T(2),T(9),0.9*mwt(i),Wgap3,Wc,Lc); % heat transfer coefficient in bottom water channel 
hc_pv2amb=1.42*(1/(0.0005/0.35+0.00005/0.85+0.0001/0.2+0.001/205+1/hwind));
%hc_filter=h_watergap(T(2),T(3),mwt(i),Wgap1,Wc,Lc); % heat transfer coefficient in the filter fluid

hc_filter=1/(Wgap2/2/kfl);

hr_pv2sky=h_radiation(Empv(1),Emss,T(8),Tss); % ideal radiation from pv to sky

%hc_pv2amb=1/(0.001/205+1/hwind);% heat conduct from pv to ambient: through the back convection

%Equations:

%Equation 1: Top channel glass
%Ac*hr_gc2sky*(Tsky-Tgc1)+Ac*hwind*(Ta-Tgc1)+Ac*(hc_wt)*(0.5*Twtin+0.5Twtout-Tgc1)+Ac*Abgc*G=0
%Ac*hr_gc2ss*(Tss-Tgc1)+Ac*hwind*(Ta-Tgc1)+Ac*(hc_wt1)*(0.5*Twtin+0.5Twtout-Tgc1)+Ac*Abgc*G=0

%Asol(1,1)=-Ac*hr_gc2sky-Ac*hwind-Ac*(hc_wt1);
Asol(1,1)=-Ac*hr_gc2ss-Ac*(hc_wt1)-Ac*hwind;
Asol(1,2)=Ac*(hc_wt1)*0.5;
Asol(1,3)=Ac*(hc_wt1)*0.5;

%Bsol(1)=-Ac*hr_gc2sky*Tsky-Ac*hwind*Ta-Ac*G1abs;
Bsol(1)=-Ac*hr_gc2ss*Tss-Ac*hwind*Ta-Ac*G1abs;
    
%Equation 2: Top water channel
%Ac*hc_wt*(Tgc1-0.5*Twtin-0.5Twtout)+Ac*hc_wt*(Tflg1-0.5*Twtin-0.5Twtout)+Ac*G2abs-mwt*cpwt*(Twtout-Twtin)=0

Asol(2,1)=Ac*hc_wt1;
Asol(2,2)=-Ac*hc_wt1*0.5-Ac*hc_wt1*0.5+0.1*mwt(i)*cpwt;
Asol(2,3)=-Ac*hc_wt1*0.5-Ac*hc_wt1*0.5-0.1*mwt(i)*cpwt;
Asol(2,4)=Ac*hc_wt1;

Bsol(2)=-Ac*G2abs;

%Equation 3: top glass of the filter
%Ac*(hc_wt)*(0.5*Twtin+0.5Twtout-Tflg1)+Ac*hc_filter*(Tflf-Tflg1)+Ac*G2abs=0%%filter


    Asol(3,2)=0.5*Ac*(hc_wt1);
    Asol(3,3)=0.5*Ac*(hc_wt1);
    Asol(3,4)=-Ac*hc_wt1-Ac*hc_filter;
    Asol(3,5)=Ac*hc_filter;
    Bsol(3)=-Ac*G2abs;
    
%Equation 4: Filter Layer
%Ac*hc_filter*(Tflg1-Tflf)+Ac*hc_filter*(Tflg2-Tflf)+Ac*G4abs-mflf*cpflf*(Tflf-Tflf_in)=0%delete mflf*cpflf*(Tflf-Tflf_in)
Asol(4,4)=Ac*hc_filter;
Asol(4,5)=-Ac*hc_filter-Ac*hc_filter;
Asol(4,6)=Ac*hc_filter;

Bsol(4)=-Ac*G4abs;
    
%Equation 5: bottom glass of the filter
%Ac*(hc_wt)*(0.5*Twtin+0.5Twtout-Tflg2)+Ac*hc_filter*(Tflf-Tflg2)+Ac*G5abs=0

    Asol(5,2)=0.5*Ac*(hc_wt2);
    Asol(5,5)=Ac*hc_filter;
    Asol(5,6)=-Ac*hc_wt2-Ac*hc_filter;
    Asol(5,9)=0.5*Ac*(hc_wt2);
    Bsol(5)=-Ac*G5abs;
    
%Equation 6: Bottom water channel
%Ac*hc_wt*(Tgc2-0.5*Twtin-0.5Twtout)+Ac*hc_wt*(Tflg2-0.5*Twtin-0.5Twtout)+Ac*G6abs-mwt*cpwt*(Twtout-Twtin)=0
Asol(6,2)=-0.5*Ac*hc_wt2-0.5*Ac*hc_wt2+0.9*mwt(i)*cpwt;
Asol(6,6)=Ac*hc_wt2;
Asol(6,7)=Ac*hc_wt2;
Asol(6,9)=-Ac*hc_wt2*0.5-Ac*hc_wt2*0.5-0.9*mwt(i)*cpwt;

Bsol(6)=-Ac*G6abs;

%Equation 7: Bottom channel glass
%Ac*hc_wt*(0.5*Twtin+0.5Twtout-Tgc2)+Ac*(hr_gc2pv+hc_gc2pv)*(Tpv-Tgc2)+Ac*G7abs+%Ac*hwind*(Ta-Tgc2)+%Ac*hr_gc2sky*(Tsky-Tgc2)=0
%Ac*hc_wt*(0.5*Twtin+0.5Twtout-Tgc2)+Ac*hr_gc2sky*(Tsky-gc2)+Asc*(hr_gc2pv+hc_gc2pv)*(Tpv-Tgc2)+Ac*G7abs=0

    Asol(7,2)=Ac*hc_wt2*0.5;
    Asol(7,7)=-Ac*hc_wt2-Ac*hr_gc22ss-Asc*(hr_gc2pv+hc_gc2pv);
    Asol(7,8)=Asc*(hr_gc2pv+hc_gc2pv);
    Asol(7,9)=Ac*hc_wt2*0.5;
    
    Bsol(7)=-Ac*G7abs-Ac*hr_gc22ss*Tss;
% %The filter far covered on the PV
% if strcmp(Distance,'far')==1
%     Asol(7,2)=Ac*hc_wt2*0.5;
%     Asol(7,7)=-Ac*hc_wt2-Ac*(hr_gc2pv+hc_gc2pv)-Ac*hwind-Ac*hr_gc2sky;
%     Asol(7,8)=Ac*(hr_gc2pv+hc_gc2pv);%delete
%     Asol(7,9)=Ac*hc_wt2*0.5;
%     
%     Bsol(7)=-Ac*G7abs-Ac*hwind*Ta-Ac*hr_gc2sky*Tsky;
% end
% 
% %The filter closely covered on the PV
% if strcmp(Distance,'close')==1
%     Asol(7,2)=Ac*hc_wt2*0.5;
%     Asol(7,7)=-Ac*hc_wt2-Ac*(hc_gc2pv+hr_gc2pv);
%     Asol(7,8)=Ac*(hc_gc2pv+hr_gc2pv);
%     Asol(7,9)=Ac*hc_wt2*0.5;
%     
%     Bsol(7)=-Ac*G7abs;%-Ac*hwind*Ta-Ac*hr_gc2sky.*Tsky;
% end

%Equation 8: PV
%Asc*(hr_gc2pv+hc_gc2pv)*(Tgc2-Tpv)+Asc*Gpvabs+Asc*hr_gc2sky*(Tsky-gc2)+Asc*hc_pv2amb*(Ta-Tpv)=0 

effm=Solarcell(1)*(1-Solarcell(6)*(T(8)-273.15-25));%module efficiency
Gpvabs=Solarcell(5)-effm*G;%waste heat in PV

Asol(8,7)=Asc*(hr_gc2pv+hc_gc2pv)-Asc*hr_gc22ss;
Asol(8,8)=-Asc*(hr_gc2pv+hc_gc2pv)-Asc*hc_pv2amb;

Bsol(8)=-Asc*Gpvabs-Asc*hc_pv2amb*Ta-Asc*hr_gc22ss*Tss;

%if strcmp(Distance,'far')==1
%     Asol(8,7)=Asc*(hc_gc2pv+hr_gc2pv);
%     Asol(8,8)=-Asc*(hc_gc2pv+hr_gc2pv)-Asc*hwind-Asc*hr_pv2sky-Asc*h_cp;
%     
%     Bsol(8)=-Asc*Gpvabs-Asc*hwind*Ta-Asc*hr_pv2sky*Tsky-Asc*h_cp*T_cp;
% end
% 
% %The filter closely covered on the PV
% %if strcmp(Distance,'close')==1
%     Asol(8,7)=Asc*(hc_gc2pv+hr_gc2pv);
%     Asol(8,8)=-Asc*(hc_gc2pv+hr_gc2pv)-Asc*h_cp;
%     
%     Bsol(8)=-Asc*Gpvabs-Asc*h_cp*T_cp;
% end



%Equation 9
% Twt=Twt
Asol(9,2)=1;
Bsol(9)=Twt;


%****LINEAR SOLVER****
X=linsolve(Asol(1:9,1:9),Bsol(1:9,1:1));
T=X.';

N=N+1;%count the iteration times
 
%end
end
%efficiencies
%Tr(i)=(0.5*T(3)+0.5*T(4)-Ta)/G;
%effthfl(i)=mwt(i)*cpwt*(0.5*T(3)+0.5*T(9)-T(2))/G/Ac;
T_ave(i)=(0.1*mwt(i)*T(3)+0.9*mwt(i)*T(9))/mwt(i);
effthfl(i)=mwt(i)*cpwt*(T_ave(i)-T(2))/G/Ac;

%effthfl(i)=mwt(i)*cpwt*(T(9)-T(2))/G/Ac;

effel(i)=effm;

%temperatures
Tsspvt(i,1)=T(1);
Tsspvt(i,2)=T(2);
Tsspvt(i,3)=T(3);
Tsspvt(i,4)=T(4);
Tsspvt(i,5)=T(5);
Tsspvt(i,6)=T(6);
Tsspvt(i,7)=T(7);
Tsspvt(i,8)=T(8);
Tsspvt(i,9)=T(9);

%energy balance

% if strcmp(Distance,'far')==1
Energy(i,1)=hr_gc2ss*(T(1)-Tsky)/G;
Energy(i,2)=hwind*(T(1)-Ta)/G;
% Energy(i,3)=hr_pv2sky*(T(8)-Tsky)/G;%far covered
Energy(i,3)=hr_gc2pv*(T(8)-T(7))/G;%close covered
% Energy(i,4)=hwind*(T(8)-Ta)/G;%Far covered
Energy(i,4)=hc_gc2pv*(T(8)-T(7))/G;%close covered
Energy(i,5)=hc_pv2amb*(T(8)-Ta)/G;%%%%%%%%%%%%%%%%CHECKKKKKKKKK
Energy(i,8)=Gpvabs/G;
end



end

% 
% %=========================calculation for the conventional PVT============
% %The structure of the  PV module

%---------------PV------------   T1
%-------------Back plate-------- T2  


% the performance of PV cell

[Trgcf,Relgcf,Abgcf]=select_filter_M('fake');

[G1pass2,G1rel2,G1abs2]=Gfilter(G0,Trgcf,Relgcf,Abgcf);%layer1:cover glass,Specrum after the cover glass%%%%CHECKKKKKKKKKK
%Gpv=ones(2002,1);
Solarcell1=solarcell(Asc,G0*Abpv,G,cell);
%simulation initials and settings
Asol1=zeros(2,2); % 2 Equations' left-hand coefficients
Bsol1=zeros(2,1); % 2 equations' right-hand coefficients
T1=zeros(1,2)+300;% T1...T2:Please see the structure of the S-PVT
X1=zeros(2,1)+300;% solution of each iteration
Ttemp1=9999;
N1=1;
error1=1e-4;
% % Solarcell1imulation initials and settings
%  Asol1=zeros(3,3); % 2 Equations' left-hand coefficients
%  Bsol1=zeros(3,1); % 2 equations' right-hand coefficients
%  Tp=zeros(1,3);%+300;% T1...T2:Please see the structure of the S-PVT
%  Xp=zeros(3,1);%+300;% solution of each iteration
%  Ttemp1=9999;
%  Np=1;
%  error1=1e-4;

 while (max(abs(X1.'-Ttemp1))>error&& N1<100)
     Ttemp1=T1;
%     
% % heat transfer coefficients calculation
hr_pv2ss=h_radiation(Empv,Emss,T1(1),Tss); % radiation from top channel to solar simulator
hwind=4*4.5+2.9*Vwind; %8.3+2.2*vwind(i-1); % 4.5+2.9*vwind, convection heat transfer due to wind


%hc_pvg2pv=1/(0.002/0.8+0.0005/0.35);% heat conduct from pv-glass to pv: through the layers of glass, eva

hc_pv2amb=1/(0.0005/0.35+0.00005/0.85+0.0001/0.2+0.001/205+1/hwind);% heat conduct from pv to abiemt: through the layers of eva, ted, glue, back plate and back convection
hc_pvg2pv=1/(0.00000004/0.8+0.000005/0.35);% heat conduct from pv-glass to pv: through the layers of glass, eva

% 
% %Equations:
% 
% %Equation 1: PV
 %Asc*hr_pv2ss*(Tss-Tpv)+Ac*hwind*(Ta-Tpv)+Asc*hc_pv2amb*(Ta-Tpv)+Asc*Gpvabs1=0
%  Asol1(1,1)=-Solarcell1(1)*Solarcell1(6);
%  Asol1(1,2)=-1;
%  Asol1(2,2)=-G;
%  Asol1(2,3)=-1;
%  Asol1(3,1)=-Asc*hr_pv2ss(2)-Ac*hwind-Asc*hc_pv2amb;
%  Asol1(3,3)=Asc;
%  Bsol1(1)=-Solarcell1(1)-Solarcell1(1)*Solarcell1(6)*298.15;
%  Bsol1(2)=-Solarcell1(5);
%  Bsol1(3)=-Asc*hr_pv2ss(2)*Tss-Asc*hwind*Ta-Asc*hc_pv2amb*Ta;
 
%  Xp=linsolve(Asol1(1:3,1:3),Bsol1(1:3,1:1));
%  Tp=Xp.';
%  Np=N1+1;%count the iteration times
%  effm1=Solarcell1(1)*(1-Solarcell1(6)*(Tp(1)-273.15-25));%module efficiency
%  Gpvabs1=Solarcell1(5)-effm1*G;%waste heat in PV
 %Equation 1: cover glass
%Ac*hr_g2sky*(Tsky-Tgpv)+Ac*hwind*(Ta-Tgpv)+Ac*hc_pvg2pv*(Tpv-Tgpv)+Ac*Abgpv*G=0

Asol1(1,1)=-Asc*hr_pv2ss(2)-Asc*hwind-Asc*hc_pvg2pv;
Asol1(1,2)=Asc*hc_pvg2pv;

Bsol1(1)=-Asc*hr_pv2ss(2)*Tss-Asc*hwind*Ta-Asc*G1abs2;


%Equation 2: PV
%Ac*hc_pvg2pv*(Tgpv-Tpv)+Ac*hc_pv2amb*(Ta-Tpv)+Ac*Gpvabs1=0
effm1=Solarcell1(1)*(1-Solarcell1(6)*(T1(2)-273.15-25));%module efficiency
Gpvabs1=Solarcell1(5)-effm1*G;%waste heat in PV


Asol1(2,1)=Asc*hc_pvg2pv;
Asol1(2,2)=-Asc*hc_pvg2pv-Asc*hc_pv2amb;

Bsol1(2)=-Asc*Gpvabs1-Asc*hc_pv2amb*Ta;



%****LINEAR SOLVER****
X1=linsolve(Asol1(1:2,1:2),Bsol1(1:2,1:1));
T1=X1.';

N1=N1+1;%count the iteration times
    
end

effel1=effm1;

%  Asol1(1,1)=-Asc*hr_pv2ss(2)-Ac*hwind-Asc*hc_pv2amb;
%  %Asol1(1,2)=Asc*hc_pv2cp;
% % 
%  Bsol1(1)=-Asc*hr_pv2ss(2)*Tss-Asc*hwind*Ta-Asc*Gpvabs1-Asc*hc_pv2amb*Ta;
% % 
% % %Equation 2: Copper Block
%  Asol(2,2)=1;
%  Bsol(9)=T_cp;
% 
% %****LINEAR SOLVER****
%  X1=linsolve(Asol1(1:3,1:3),Bsol1(1:3,1:1));
%  Tp=Xp.';
% 
 N1=N1+1;%count the iteration times
%     
% end
% effel1=effm1;
% 
 %temperatures
Tcpvt(1)=T1(1);
Tcpvt(2)=T1(2);
% 
% 
% %energy balance
 Energy1(1)=hr_pv2ss(2)*(T1(1)-Tsky)/G;
 Energy1(2)=hwind*(T1(1)-Ta)/G;
 Energy1(3)=hc_pv2amb*(T1(2)-Ta)/G;
 Energy1(4)=Gpvabs1/G;
% 
% 
% 
%calculate the special ponit value:the electrical, thermal improvement at 65C outlet water
% mfi=interp1(Tsspvt(:,4),mwt',338,'linear')/Ac;
% 
% P(1)=interp1(mwt',Tsspvt(:,7),mfi*Ac,'linear')-273.15;
% 
% P(2)=Tcpvt(2)-273.15;
% 
% P(3)=P(2)-P(1);
% 
% P(4)=interp1(mwt,effel,mfi*Ac,'linear');
% 
% P(5)=effel1;
% 
% P(6)=P(4)-P(5);

 