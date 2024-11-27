%Anil Cengiz - Oct 19,2024

clear all
clc
close all


%% load PKPD parameters
load('data_parameters.mat')

ode_options = odeset('MaxStep',24,'RelTol',1e-8,'AbsTol',1e-10);

%%
%Initial Conditions for Body Weight Effects
bw_pld_ic=0;
bw_pln_ic=0;
bw_s_ic=0;


%Semaglutide
dosing_interval=[]*24;

max_time=24*7*52*10;


tau=7*24;
[t_semag,conc]=conc_long(semag,tau,max_time);
c_semag=1e6*conc'/semag.molar_mass;
c_semag_av=filter(ones(1,7*24)/(7*24),1,c_semag);
c_av_ss=c_semag_av(end)*ones(length(c_semag_av),1);


[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag,c_semag,semag),t_semag,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
semag_bw=semag.bw0+y(:,1)+y(:,2)-y(:,3)-semag.bw0*semag.emaxi*c_semag./(semag.ec50i+c_semag);
semag_bw=100*(semag_bw-semag.bw0)/semag.bw0;


[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag,c_semag_av,semag),t_semag,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
semag_bw_av=semag.bw0+y(:,1)+y(:,2)-y(:,3)-semag.bw0*semag.emaxi*c_semag_av./(semag.ec50i+c_semag_av);
semag_bw_av=100*(semag_bw_av-semag.bw0)/semag.bw0;


[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag,c_semag,semag),t_semag,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
semag_bw_av_ss=semag.bw0+y(:,1)+y(:,2)-y(:,3)-semag.bw0*semag.emaxi*c_av_ss./(semag.ec50i+c_av_ss);
semag_bw_av_ss=100*(semag_bw_av_ss-semag.bw0)/semag.bw0;


t_semag=t_semag/(7*24);





col_vec=["#e41a1c","#377eb8","#4daf4a","#984ea3"];


%% Graphing

linwid=2;
hfig=figure()
hold on;

% thin=10;
% t_semag=t_semag(1:thin:end);
% semag_bw=semag_bw(1:thin:end);
% semag_bw_av=semag_bw_av(1:thin:end);
plot(t_semag,semag_bw,'Linewidth',linwid, 'Color',col_vec(1));
plot(t_semag,semag_bw_av,'--','Linewidth',linwid, 'Color',col_vec(2));
% plot(t_semag,semag_bw_av_ss,'Linewidth',linwid, 'Color',col_vec(3));


%plot(dosing_interval/24,semag_ss*100,'LineWidth',linwid,'Color',col_vec(1));

%plot(dosing_interval/24,tirz_5*100,'LineWidth',linwid,'Color',col_vec(2));

%plot(dosing_interval/24,tirz_10*100,'LineWidth',linwid,'Color',col_vec(3));

%plot(dosing_interval/24,tirz_15*100,'LineWidth',linwid,'Color',col_vec(4));
%legend("Semaglutide 2.4 mg", "Tirzepatide 5mg", "Tirzepatide 10mg", "Tirzepatide 15mg", "Location","Southeast")
% legend("$C(t)$","$C_{av}(t)$","$C_{av,ss}$")
legend("$C$ used in PD model","$C_{\textup{avg}}$ used in PD model")
xlabel("Time (Weeks)")
ylabel("Percent Change in Body Weight")
grid on
grid minor
xlim([0,80])
%ylabel("$\%\Delta$ in BW ($\%$)")


%fname="smooth_no_hill9_ss_comp";
%fname="fmincon_ss_comp_better_initial_4";
%fname="psort_ss_comp_better_initial_9555";
fname="figSemagCheck"

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
exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(hfig,strcat(fname,'.fig'))

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


