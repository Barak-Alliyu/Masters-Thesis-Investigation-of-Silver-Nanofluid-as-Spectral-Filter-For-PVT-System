%Figure1: spectrum distributions
figure;
x0=0;
y0=0;
width=560*1.5;
height=420*1.5;
set(gcf,'position',[x0,y0,width,height]);

x_G=xlsread('data', 'ambient conditions', 'B2:B2003');
SR_CdTe=xlsread('data', 'solar cells', 'D3:D2004');
SR_Sil=xlsread('data', 'solar cells', 'I3:I2004');

Plotspectr(1)=plot(x_G,G0,'-','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{Original Spectrum}');hold on;
Plotspectr(2)=plot(x_G,G4pass,'-','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{Filtered Spectrum}');hold on;
Plotspectr(3)=plot(x_G,SR_CdTe,'--','Color',[1 0 0],'linewidth',2,'DisplayName','{CdTe Spectral Response}');
%Plotspectr(3)=plot(x_G,SR_Sil,'--','Color',[0 0 0],'linewidth',2,'DisplayName','{Silicon Spectral Response}');

axis([0 2500 0 1.8]);
xlabel('{Wavelength} [nm]','FontSize',16);
ylabel('{Spectral intensity [W/m2/nm]   Spectral Response [W/A]}','FontSize',14);
legenspectr=legend(Plotspectr,'Location','northeast');legend boxon;
set(legenspectr,'FontSize',18);
set(gca,'Layer','Top','XColor','k','YColor','k','linewidth',0.5,'FontSize',16,'XTick',[0:500:2500],'YTick',[0:0.2:1.8],'TickDir','in','box','off');
set(gcf,'color','white')
box on;



% % %figure 2-1: Temperature distributions via mass flow rate

figure;
x0=0;
y0=800;
width=560*1.5;
height=420*1.5;
set(gcf,'position',[x0,y0,width,height]);


Plottem(1)=plot(mfl/Ac,Tsspvt(:,7)-273,'-o','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{Tpv-SS}');hold on;
Plottem(2)=plot(mfl/Ac,Tsspvt(:,4)-273,'--o','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{Tout-SS}');hold on;

Tcpvtm=ones(Nm(2),1)*Tcpvt(2);
Plottem(3)=plot(mfl/Ac,Tcpvtm-273,'-o','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{Tpv-C}');hold on;
%  Plottem(5)=plot(mc,Tcpvt(:,5)-273,'--o','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{Tout-C}');hold on;



grid on;

axis([0 0.003 27 110]);
xlabel('{M}_{fluid} [-]','FontSize',16);
ylabel('{Temperature [C]}','FontSize',16);
legendtem=legend(Plottem,'Location','northeast');legend boxon;
set(legendtem,'FontSize',12);
set(gca,'Layer','Top','XColor','k','YColor','k','linewidth',0.5,'FontSize',16,'XTick',[0:0.0005:0.006],'YTick',[30:10:200],'TickDir','in','box','off');
set(gcf,'color','white')
box on;



%figure 3-1 efficiencies via mass flow rate

figure;
x0=500;
y0=0;
width=560*1.5;
height=420*1.5;
set(gcf,'position',[x0,y0,width,height]);


Ploteff1(1)=plot(mfl/Ac,effthfl,'--o','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{\it\eta}_{t-Spectral splitting PVT}');hold on;
%Ploteff1(2)=plot(mfl,effth1,'--o','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{\it\eta}_{t-Conventional PVT}');hold on;
%Ploteff1(3)=plot(mc,effthfl,'-','Color',[0 0 0],'linewidth',2,'DisplayName','{\it\eta}_{t-filter}');hold on;
%Plotspectr(3)=plot(Tr,effthcool,'-o','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{\it\eta}_{t-coolant}');hold on;
Ploteff1(2)=plot(mfl/Ac,effel,'-o','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{\it\eta}_{el-Spectral splitting PVT}');

effel1m=ones(Nm(2),1)*effel1;
Ploteff1(3)=plot(mfl/Ac,effel1m,'-o','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{\it\eta}_{el-Conventional PV}');
grid on;

axis([0 0.003 0 0.5]);
xlabel('{M}_{fluid} [-]','FontSize',16);
ylabel('{Efficiency [-]}','FontSize',16);
legendeff1=legend(Ploteff1,'Location','northeast');legend boxon;
set(legendeff1,'FontSize',12);
set(gca,'Layer','Top','XColor','k','YColor','k','linewidth',0.5,'FontSize',16,'XTick',[0:0.0005:0.006],'YTick',[0:0.05:1],'TickDir','in','box','off');
set(gcf,'color','white')
box on;



%figure 3: energy balance of the PV

figure;
x0=1000;
y0=0;
width=560*1.5;
height=420*1.5;
set(gcf,'position',[x0,y0,width,height]);


Plotenergy(1)=plot(mfl/Ac,Energy(:,8),'-','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{Q}_{absorb-SS}');hold on;
Plotenergy(2)=plot(mfl/Ac,Energy(:,3),'--','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{Q}_{loss-radiation-top-SS}');hold on;
Plotenergy(3)=plot(mfl/Ac,Energy(:,4),'--p','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{Q}_{loss-covection-top-SS}');hold on;
Plotenergy(4)=plot(mfl/Ac,Energy(:,5),'--o','Color',[0.8500 0.3250 0.0980],'linewidth',2,'DisplayName','{Q}_{loss-covection-back-SS}');hold on;



Plotenergy(5)=plot(mfl/Ac,ones(Nm(2),1)*Energy1(4),'-','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{Q}_{absorb-C}');hold on;
Plotenergy(6)=plot(mfl/Ac,ones(Nm(2),1)*Energy1(1),'--','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{Q}_{loss-radiation-top-C}');hold on;
Plotenergy(7)=plot(mfl/Ac,ones(Nm(2),1)*Energy1(2),'--p','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{Q}_{loss-covection-top-C}');hold on;
Plotenergy(8)=plot(mfl/Ac,ones(Nm(2),1)*Energy1(3),'--o','Color',[0 0.4470 0.7410],'linewidth',2,'DisplayName','{Q}_{loss-covection-back-C}');


grid on;

axis([0.001 0.003 0 1]);
xlabel('{M}_{fluid} [-]','FontSize',16);
ylabel('{Energy flow [-]}','FontSize',16);
legendenergy=legend(Plotenergy,'Location','southwest');legend boxon;
set(legendenergy,'FontSize',12);
set(gca,'Layer','Top','XColor','k','YColor','k','linewidth',0.5,'FontSize',16,'XTick',[0:0.0005:0.006],'YTick',[-1:0.1:1],'TickDir','in','box','off');
set(gcf,'color','white')
box on;


%figure 4: energy flow of the solar energy
figure;
G_abs_filter=(G2abs+G3abs+G4abs)/G;
G_abs_pv=Gpvabs/G;
G_el=effm;
G_rel=1-G_abs_filter-G_abs_pv-G_el;
Gflow=[G_abs_filter,G_abs_pv,G_el,G_rel];
labels={'Filter','PV-th','PV-el','Reflective Loss'};

p=pie(Gflow);

legend(labels,'Location','southoutside','Orientation','horizontal');
set(gcf,'color','white')

title('Solar Energy Flow','FontSize',15);



