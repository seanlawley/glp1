%Anil Cengiz - Nov 15, 2024
%This scripts simulates the fitted semaglutide model and generates
% %timecourses
%
close all
clear all
% clc

%% load PKPD parameters
load('data_parameters.mat')

%% load data from Wilding et al 2021

W1=  readmatrix("wilding_2021_data.csv");

%Clean Up of Data File
W1=W1(3:end,:);

W=W1;
W(:,[1,2,5,6,7,8,11,12])=[[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];W(1:end-1,[1,2,5,6,7,8,11,12])];

%Assign error bars and data points
t_pl_top=W(:,1);
pl_top=W(:,2);
t_pl_data=W(:,3);
pl_data=W(:,4);

t_pl_bot=W(:,5);
pl_bot=W(:,6);

t_semag_top=W(:,7);
semag_top=W(:,8);

t_semag_data=W(:,9);
semag_data=W(:,10);

t_semag_bot=W(:,11);
semag_bot=W(:,12);

%%

max_time=7*24*74;
tau=7*24;
tspan=0:max_time;
bw0=104.8;

% Semaglutide
[t_semag,c_semag]=conc_time_courses_mod(semag,2,7*24,0,80*7*24);
%c_semag=c_semag(2,:)'*1e6/par.molar_mass;
t_semag=t_semag';
%Placebo Semaglutide
t_pl=t_semag;
c_pl=0*t_pl';

y0=[0,0,0];
ode_options = odeset('MaxStep',24);
%%
%Computes the slow placebo effect, fast placebo effect, and slow drug
%effect
[~,y_pl_ode] = ode45(@(t,y)odefcn(t,y,semag,t_pl,c_pl,bw0),t_pl,y0,ode_options);
[~,y_semag_ode] = ode45(@(t,y)odefcn(t,y,semag,t_semag,c_semag,bw0),t_semag,y0,ode_options);
%%
%Computes time courses for drug effects
%y=bw0+y_ode(:,1) + y_ode(:,2) - y_ode(:,3)-alg_func
%y_ode(:,1):fast placebo effect
%y_ode(:,2):slow placebo effect
%y_ode(:,3):slow drug effect
%alg_func:fast drug effect

emaxi=semag.emaxi;
n1=semag.n1;
ec50i=semag.ec50i;
y_pl=bw0+y_pl_ode(:,1)+y_pl_ode(:,2);
y_semag=bw0+y_semag_ode(:,1)+y_semag_ode(:,2)-y_semag_ode(:,3)-bw0*emaxi.*(c_semag.^n1)./(ec50i+(c_semag.^n1));

%Convert body weights into change in body weight (%)
y_pl=100*(y_pl-bw0)/bw0;
y_semag=100*(y_semag-bw0)/bw0;

%%
%Plot time course figures
linwid=1.5;
hfig=figure()
hold on


e1=errorbar(t_pl_data,pl_data,pl_data-pl_bot,pl_top-pl_data,'s','MarkerSize',10,'MarkerEdgeColor','black')
e1.Color='black';
e1.LineWidth=linwid;

e2=errorbar(t_semag_data,semag_data,semag_data-semag_bot,semag_top-semag_data,'^','MarkerSize',10,'MarkerEdgeColor','#e41a1c')
e2.Color='#e41a1c';
e2.LineWidth=linwid;

plot(t_pl/(7*24),y_pl,'Color','black','LineWidth',linwid)
plot(t_semag/(7*24),y_semag,'Color','#e41a1c','LineWidth',linwid)

% title("Semaglutide Efficacy with 16 Week Dose Escalation Period")
% yline(0,"LineStyle","-.")
% text(55,0.5,"Baseline Body Weight=104.8 kg")
legend("Placebo","Semaglutide 2.4 mg",'Location','southwest')
grid on
grid minor

xlabel("Time (Weeks)")
ylabel("Percent Change in Body Weight")
ylim([-25,1])
xlim([0,70])


fname=sprintf("figSemagFit");

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

function dydt=odefcn(t,y,semag,t_data,c_data,bw0)
%y(1):fast placebo effect
%y(2):slow placebo effect
%y(3):slow drug effect

%ODE Solving function
dydt=zeros(3,1);
%Load PD parameters
kout=semag.kout;
emaxi=semag.emaxi;
koutpl=semag.koutpl;
plbmax=semag.plbmax;
tfac=semag.tfac;
emaxs=semag.emaxs;
ec50i=semag.ec50i;
ec50s=semag.ec50s;

newss=semag.newss;
kout_s_plcb=semag.kout_s_plcb;

n1=semag.n1;
n2=semag.n2;

%Interpolate drug concentration
c=interp1(t_data,c_data,t);
%Solve ODEs
dydt(1)=-koutpl*(y(1)+bw0*plbmax*exp(-tfac*t));
dydt(2)=-kout_s_plcb*(y(2)-bw0*(newss-1));
dydt(3)=-kout*(y(3)-newss*bw0*emaxs*((c^n2)/(ec50s+(c^n2))));
end
