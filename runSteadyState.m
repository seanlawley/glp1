%Anil Cengiz - Oct 19,2024

%This script runs parallel loops to run long time course simulations and evaluate the steady-state numerical estimation of body weight

clear all
clc
close all
tic

%% load PKPD parameters
load('data_parameters.mat')


col_vec=["#e41a1c","#377eb8","#4daf4a","#984ea3"];
ode_options = odeset('MaxStep',24,'RelTol',1e-8,'AbsTol',1e-10);

%%
%Initial Conditions for Body Weight Effects 
bw_pld_ic=0;
bw_pln_ic=0;
bw_s_ic=0;


%Semaglutide
dosing_interval=[7:28]*24;

max_time=24*7*52*10;
semag_av=zeros(size(dosing_interval));
tirz_av_5=zeros(size(dosing_interval));
tirz_av_10=zeros(size(dosing_interval));
tirz_av_15=zeros(size(dosing_interval));
tirz5=tirz;
tirz10=tirz;
tirz15=tirz;
tirz5.weekly_dose=5;
tirz10.weekly_dose=10;
tirz15.weekly_dose=15;



parfor j=1:length(dosing_interval)
    tau=dosing_interval(j);
    [t_semag,conc]=conc_long(semag,tau,max_time);
    c_semag=1e6*conc'/semag.molar_mass;

     [t_tirz_5,conc]=conc_long(tirz5,tau,max_time);
    c_tirz_5=1e6*conc'/tirz.molar_mass;

     [t_tirz_10,conc]=conc_long(tirz10,tau,max_time);
    c_tirz_10=1e6*conc'/tirz.molar_mass;

     [t_tirz_15,conc]=conc_long(tirz15,tau,max_time);
    c_tirz_15=1e6*conc'/tirz.molar_mass;

[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag,c_semag,semag),t_semag,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
semag_bw=semag.bw0+y(:,1)+y(:,2)-y(:,3)-semag.bw0*semag.emaxi*c_semag./(semag.ec50i+c_semag);
semag_bw=100*(semag_bw-semag.bw0)/semag.bw0;
semag_av(j)=trapz(t_semag,semag_bw)/t_semag(end);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz_5,c_tirz_5,tirz),t_tirz_5,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
tirz_bw_5=tirz.bw0+y(:,1)+y(:,2)-y(:,3)-tirz.bw0*tirz.emaxi*c_tirz_5./(tirz.ec50i+c_tirz_5);
tirz_bw_5=100*(tirz_bw_5-tirz.bw0)/tirz.bw0;
tirz_av_5(j)=trapz(t_tirz_5,tirz_bw_5)/t_tirz_5(end);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz_10,c_tirz_10,tirz),t_tirz_10,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
tirz_bw_10=tirz.bw0+y(:,1)+y(:,2)-y(:,3)-tirz.bw0*tirz.emaxi*c_tirz_10./(tirz.ec50i+c_tirz_10);
tirz_bw_10=100*(tirz_bw_10-tirz.bw0)/tirz.bw0;
tirz_av_10(j)=trapz(t_tirz_10,tirz_bw_10)/t_tirz_10(end);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz_15,c_tirz_15,tirz),t_tirz_15,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
tirz_bw_15=tirz.bw0+y(:,1)+y(:,2)-y(:,3)-tirz.bw0*tirz.emaxi*c_tirz_15./(tirz.ec50i+c_tirz_15);
tirz_bw_15=100*(tirz_bw_15-tirz.bw0)/tirz.bw0;
tirz_av_15(j)=trapz(t_tirz_15,tirz_bw_15)/t_tirz_15(end);


end 

semag_plc=semag.newss-1;
D=1e6*2.4/semag.molar_mass;
%semag=(semag.newss-1)-(semag.emaxi+semag.emaxs*semag.newss)*D*semag.F./(D*semag.F+semag.ec50*semag.cl*tau);
tau=dosing_interval;
semag_ss=(semag.newss-1)-semag.emaxi*D*semag.F./(D*semag.F+semag.ec50i*semag.V*semag.ke*tau)-semag.emaxs*D*semag.F./(D*semag.F+semag.ec50s*semag.V*semag.ke*tau)

D5=1e6*5/tirz.molar_mass;
D10=1e6*10/tirz.molar_mass;
D15=1e6*15/tirz.molar_mass;
tirz_plc=tirz.newss-1;
tirz_5=(tirz.newss-1)-tirz.emaxi*D5*tirz.F./(D5*tirz.F+tirz.ec50i*tirz.cl*tau)-tirz.emaxs*D5*tirz.F./(D5*tirz.F+tirz.ec50s*tirz.cl*tau);
tirz_10=(tirz.newss-1)-tirz.emaxi*D10*tirz.F./(D10*tirz.F+tirz.ec50i*tirz.cl*tau)-tirz.emaxs*D10*tirz.F./(D10*tirz.F+tirz.ec50s*tirz.cl*tau)
tirz_15=(tirz.newss-1)-tirz.emaxi*D15*tirz.F./(D15*tirz.F+tirz.ec50i*tirz.cl*tau)-tirz.emaxs*D15*tirz.F./(D15*tirz.F+tirz.ec50s*tirz.cl*tau)


toc

%% Graphing

linwid=1.5;
hfig=figure()
hold on;
scatter(dosing_interval/(24),semag_av,80,'MarkerEdgeColor',col_vec(1),'MarkerFaceColor',col_vec(1));
scatter(dosing_interval/(24),tirz_av_5,80,'MarkerEdgeColor',col_vec(2),'MarkerFaceColor',col_vec(2));
scatter(dosing_interval/(24),tirz_av_10,80,'MarkerEdgeColor',col_vec(3),'MarkerFaceColor',col_vec(3));
scatter(dosing_interval/(24),tirz_av_15,80,'MarkerEdgeColor',col_vec(4),'MarkerFaceColor',col_vec(4));


%plot(dosing_interval/24,semag_ss*100,'LineWidth',linwid,'Color',col_vec(1));

%plot(dosing_interval/24,tirz_5*100,'LineWidth',linwid,'Color',col_vec(2));

%plot(dosing_interval/24,tirz_10*100,'LineWidth',linwid,'Color',col_vec(3));

%plot(dosing_interval/24,tirz_15*100,'LineWidth',linwid,'Color',col_vec(4));
legend("Semaglutide 2.4 mg", "Tirzepatide 5mg", "Tirzepatide 10mg", "Tirzepatide 15mg", "Location","Southeast")

xlabel("Dosing Interval (Days)")
ylabel("\% Change in BW at SS")
grid on 
grid minor
% xlim([4,29])
%ylabel("$\%\Delta$ in BW ($\%$)")


%% format and save figure

fname="optimal_pset_ss"

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

legend boxoff
%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
% exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
% saveas(hfig,strcat(fname,'.fig'))

%% save data

% timestamp = char(datetime('now'));
% timestamp = strrep(timestamp,' ','-');
% timestamp = strrep(timestamp,':','');
% filename = strcat(['data_SteadyState_',timestamp,'.mat']);
filename = strcat(['data_SteadyState.mat']);
save(filename);

%%

function bw=body_weight(bw_0, bw_s, bw_i, bw_pld, bw_pln)
bw=bw_0-bw_s-bw_i+bw_pld+bw_pln;
end

function dydt = ode_solver(t,y,tc,c,par)
   bw_pld=y(1);
   bw_pln=y(2);
   bw_s=y(3);

  cavg=interp1(tc,c,t);

  kin=par.bw0*par.kout*par.newss; 
    
  dydt = zeros(3,1);
  dbw_pld=-par.koutpl*(bw_pld+par.bw0*par.plbmax*exp(-par.tfac*t));
  dbw_pln=-par.kout_s_plcb*(bw_pln-par.bw0*(par.newss-1));
  
  dbw_s=kin*(par.emaxs*cavg/(par.ec50s+cavg))-par.kout*bw_s;
  %dbw_s=kin*(0.01*emaxs*cavg/(ec50+cavg))-kout*bw_s;

  dydt(1) = dbw_pld;
  dydt(2) =dbw_pln;
  dydt(3)=dbw_s;
end


