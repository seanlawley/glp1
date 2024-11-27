%Anil Cengiz - Nov 25,2024

%This script runs parallel loops to compute efficacy at steady state for
%different dosing intervals and different emaxi emaxs values.

%It generates a plot of results superimposed and saves the results in the
%semag_ss.mat file.

clear all
clc
close all

col_vec=["#e41a1c","#377eb8","#4daf4a","#984ea3"];

p_estimate=[0.0419,0.2257];


%Semaglutide PK 
semag.name="semaglutide";
semag.ka=0.0296; %absorption rate (1/hour)
semag.V=12.1; %TABLE S2A
semag.ke=0.0488/semag.V; %TABLE S2A 
semag.F=1;
semag.weekly_dose=2.4;
semag.bw0=104.8;
semag.molar_mass=4114;

%Set PD parameters to estimated values


semag.F=1;
semag.molar_mass=4114; % g/mol 
semag.kout= 0.0319/(7*24);
semag.kout_s_plcb= 0.0319/(7*24);

semag.emaxi=0.01*3.99;
semag.koutpl= 0.00885/(7*24);
semag.plbmax=0.01*34.4;
semag.tfac=0.0794/(7*24);
semag.emaxs=0.01*26.2;
semag.ec50i=48.0;
semag.ec50s=48.0;

semag.newss=0.990;
semag.bw0=104.8;
semag.weekly_dose=2.4;

emaxi_vec=p_estimate(1)*[0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3,1.4, 1.5];
emaxs_vec=p_estimate(2)*[0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3,1.4, 1.5];

ode_options = odeset('MaxStep',24,'RelTol',1e-8,'AbsTol',1e-10);
%%
%Initial Conditions for Body Weight Effects 
bw_pld_ic=0;
bw_pln_ic=0;
bw_s_ic=0;

%Semaglutide
dosing_interval=[7:28]*24;
dosing_interval=[5:1:28]*24;

holder=zeros(length(emaxi_vec),length(emaxs_vec),length(dosing_interval));

max_time=24*7*52*10;

linwid=1.5;
hfig=figure()
hold on;

for i=1:length(emaxi_vec)
 for  k=1:length(emaxs_vec)
   semag.emaxi=emaxi_vec(i);
   semag.emaxs=emaxs_vec(k);

semag_av=zeros(size(dosing_interval));

parfor j=1:length(dosing_interval)
    tau=dosing_interval(j);
    [t_semag,conc]=conc_long(semag,tau,max_time);
    c_semag=1e6*conc'/semag.molar_mass;


[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag,c_semag,semag),t_semag,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
semag_bw=semag.bw0+y(:,1)+y(:,2)-y(:,3)-semag.bw0*semag.emaxi*c_semag./(semag.ec50i+c_semag);
semag_bw=100*(semag_bw-semag.bw0)/semag.bw0;
semag_av(j)=trapz(t_semag,semag_bw)/t_semag(end);

end 

holder(i,k,:)=semag_av;


%% Graphing


scatter(dosing_interval/(24),semag_av,80,'MarkerEdgeColor',col_vec(1),'MarkerFaceColor',col_vec(1),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
 end
end
legend("Semaglutide 2.4 mg", "Location","Southeast")
xlabel("Dosing Interval (Days)")
ylabel("\% Change in BW at SS")
grid on 
grid minor
xlim([4,29])

fname="ss_semag_bw_dosing_ints"

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


