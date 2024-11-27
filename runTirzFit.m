%Anil Cengiz - Oct 15, 2024
%This scripts simulates the fitted tirzepatide model and generates
% %timecourses
% 
close all
clear all

%% load PKPD parameters
load('data_parameters.mat')

%% load data from Jastreboff et al 2022

JM=  readmatrix("jastreboff_2022_data_mod.csv");
JM=JM(3:end-1,:);
pl_top=JM(:,4);
pl_bot=JM(:,6);

t_5_top=JM(:,10);
t_5_bot=JM(:,12);

t_10_top=JM(:,16);
t_10_bot=JM(:,18);

t_15_top=JM(:,22);
t_15_bot=JM(:,24);

max_time=7*24*74;
tau=7*24;
tspan=0:max_time;
bw0=104.8;

%5mg Tirzepatide
[t_5,c_5]=conc_time_courses(tirz,2,7*24,5);
c_5=c_5(2,:)'*1e6/tirz.molar_mass;
%10mg Tirzepatide
[t_10,c_10]=conc_time_courses(tirz,2,7*24,10);
c_10=c_10(2,:)'*1e6/tirz.molar_mass;
%15mg Tirzepatide
[t_15,c_15]=conc_time_courses(tirz,2,7*24,15);
c_15=c_15(2,:)'*1e6/tirz.molar_mass;
%Placebo Tirzpatide
t_pl=t_5; 
c_pl=0*t_pl';

y0=[0,0,0];
ode_options = odeset('MaxStep',24);
%%
%Computes the slow placebo effect, fast placebo effect, and slow drug
%effect
[~,y_pl_ode] = ode45(@(t,y)odefcn(t,y,tirz,t_pl,c_pl,bw0),t_pl,y0,ode_options);
[~,y_5_ode] = ode45(@(t,y)odefcn(t,y,tirz,t_5,c_5,bw0),t_5,y0,ode_options);
[~,y_10_ode] = ode45(@(t,y)odefcn(t,y,tirz,t_10,c_10,bw0),t_10,y0,ode_options);
[~,y_15_ode] = ode45(@(t,y)odefcn(t,y,tirz,t_15,c_15,bw0),t_15,y0,ode_options);

emaxi=tirz.emaxi;
n1=tirz.n1;
ec50i=tirz.ec50i;
y_pl=bw0+y_pl_ode(:,1)+y_pl_ode(:,2);
y_5=bw0+y_5_ode(:,1)+y_5_ode(:,2)-y_5_ode(:,3)-bw0*emaxi.*(c_5.^n1)./(ec50i+(c_5.^n1));
y_10=bw0+y_10_ode(:,1)+y_10_ode(:,2)-y_10_ode(:,3)-bw0*emaxi*(c_10.^n1)./(ec50i+(c_10.^n1));
y_15=bw0+y_15_ode(:,1)+y_15_ode(:,2)-y_15_ode(:,3)-bw0*emaxi*(c_15.^n1)./(ec50i+(c_15.^n1));


%Convert body weights into change in body weight (%)
y_pl=100*(y_pl-bw0)/bw0;
y_5=100*(y_5-bw0)/bw0;
y_10=100*(y_10-bw0)/bw0;
y_15=100*(y_15-bw0)/bw0;


ind=find(t_pl/(7*24)>40,1);
% display([y_pl(ind),y_5(ind),y_10(ind),y_15(ind)]);

ind=find(t_pl/(7*24)>52,1);
% display([y_pl(ind),y_5(ind),y_10(ind),y_15(ind)]);

%% save data

% timestamp = char(datetime('now'));
% timestamp = strrep(timestamp,' ','-');
% timestamp = strrep(timestamp,':','');
% filename = strcat(['data_TirzFit_',timestamp,'.mat']);
filename = strcat(['data_TirzFit.mat']);
save(filename);

%% Graphing

%%
ms=10; % marker size
%Plot time course figures 
linwid=1.5;
hfig=figure();
hold on

e1=errorbar(JM(:,1),JM(:,2),JM(:,2)-pl_bot,pl_top-JM(:,2),'s','MarkerSize',10,'MarkerEdgeColor','black');
e1.Color='black';
e1.LineWidth=linwid;
e2=errorbar(JM(:,7),JM(:,8),JM(:,8)-t_5_bot,t_5_top-JM(:,8),'o','MarkerSize',10,'MarkerEdgeColor','#377eb8');
e2.Color='#377eb8';
e2.LineWidth=linwid;
e3=errorbar(JM(:,13),JM(:,14),JM(:,14)-t_10_bot,t_10_top-JM(:,14),'v','MarkerSize',10,'MarkerEdgeColor','#4daf4a');
e3.Color='#4daf4a';
e3.LineWidth=linwid;
e4=errorbar(JM(:,19),JM(:,20),JM(:,20)-t_15_bot,t_15_top-JM(:,20),'diamond','MarkerSize',10,'MarkerEdgeColor','#984ea3');
e4.Color='#984ea3';
e4.LineWidth=linwid;

plot(t_pl/(7*24),y_pl,'Color','black','LineWidth',linwid)
plot(t_5/(7*24),y_5,'Color','#377eb8','LineWidth',linwid)
plot(t_10/(7*24),y_10,'Color','#4daf4a','LineWidth',linwid)
plot(t_15/(7*24),y_15,'Color','#984ea3','LineWidth',linwid)

legend("Placebo","Tirzepatide 5mg","Tirzepatide 10mg","Tirzepatide 15mg",'Location','southwest')
grid on 
grid minor

xlabel("Time (Weeks)")
ylabel("Percent Change in Body Weight")
ylim([-25,1])
xlim([0,80])


fname=sprintf("figTirzFit");


picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')

% text(-10.5, 0.1, 'A', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');


legend boxoff
%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(hfig,strcat(fname,'.fig'))


function dydt=odefcn(t,y,p_est,t_data,c_data,bw0)
%y(1):fast placebo effect
%y(2):slow placebo effect
%y(3):slow drug effect

%ODE Solving function
dydt=zeros(3,1);
%Load PD parameters
kout=p_est.kout;
emaxi=p_est.emaxi;
koutpl=p_est.koutpl;
plbmax=p_est.plbmax;
tfac=p_est.tfac;
emaxs=p_est.emaxs;
ec50i=p_est.ec50i;
ec50s=p_est.ec50s;

newss=p_est.newss;
kout_s_plcb=p_est.kout_s_plcb;

n1=p_est.n1;
n2=p_est.n2;

%Interpolate drug concentration
c=interp1(t_data,c_data,t);
%Solve ODEs
dydt(1)=-koutpl*(y(1)+bw0*plbmax*exp(-tfac*t));
dydt(2)=-kout_s_plcb*(y(2)-bw0*(newss-1));
dydt(3)=-kout*(y(3)-newss*bw0*emaxs*((c^n2)/(ec50s+(c^n2))));
end
