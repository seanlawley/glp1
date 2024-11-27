clear all
% clc
close all

%% load PKPD parameters
load('data_parameters.mat')

semag.weekly_dose=2.4;

ode_options = odeset('MaxStep',24,'RelTol',1e-8,'AbsTol',1e-10);
%%
%Initial Conditions for Body Weight Effects 
bw_pld_ic=0;
bw_pln_ic=0;
bw_s_ic=0;

tau=7*24;
dos_switch_time=120*7*24; %When we switch to reduced dosing
tpoint=dos_switch_time*2; %For table creation

[t_semag_pl1,semag_conc_pl1]=conc_time_courses_mod(semag,0,tau,0,dos_switch_time);
[t_tirz5_pl1,tirz5_conc_pl1]=conc_time_courses_mod(tirz,0,tau,5,dos_switch_time);
[t_tirz10_pl1,tirz10_conc_pl1]=conc_time_courses_mod(tirz,0,tau,10,dos_switch_time);
[t_tirz15_pl1,tirz15_conc_pl1]=conc_time_courses_mod(tirz,0,tau,15,dos_switch_time);

[t_semag_pl2,semag_conc_pl2]=conc_time_courses_mod(semag,1,tau,0,dos_switch_time);
[t_tirz5_pl2,tirz5_conc_pl2]=conc_time_courses_mod(tirz,1,tau,5,dos_switch_time);
[t_tirz10_pl2,tirz10_conc_pl2]=conc_time_courses_mod(tirz,1,tau,10,dos_switch_time);
[t_tirz15_pl2,tirz15_conc_pl2]=conc_time_courses_mod(tirz,1,tau,15,dos_switch_time);

[t_semag_7,semag_conc_7]=conc_time_courses_mod(semag,2,tau,0,dos_switch_time);
[t_tirz5_7,tirz5_conc_7]=conc_time_courses_mod(tirz,2,tau,5,dos_switch_time);
[t_tirz10_7,tirz10_conc_7]=conc_time_courses_mod(tirz,2,tau,10,dos_switch_time);
[t_tirz15_7,tirz15_conc_7]=conc_time_courses_mod(tirz,2,tau,15,dos_switch_time);

tau=10*24;
[t_semag_10,semag_conc_10]=conc_time_courses_mod(semag,2,tau,0,dos_switch_time);
[t_tirz5_10,tirz5_conc_10]=conc_time_courses_mod(tirz,2,tau,5,dos_switch_time);
[t_tirz10_10,tirz10_conc_10]=conc_time_courses_mod(tirz,2,tau,10,dos_switch_time);
[t_tirz15_10,tirz15_conc_10]=conc_time_courses_mod(tirz,2,tau,15,dos_switch_time);

tau=14*24;
[t_semag_14,semag_conc_14]=conc_time_courses_mod(semag,2,tau,0,dos_switch_time);
[t_tirz5_14,tirz5_conc_14]=conc_time_courses_mod(tirz,2,tau,5,dos_switch_time);
[t_tirz10_14,tirz10_conc_14]=conc_time_courses_mod(tirz,2,tau,10,dos_switch_time);
[t_tirz15_14,tirz15_conc_14]=conc_time_courses_mod(tirz,2,tau,15,dos_switch_time);

tau=21*24;
[t_semag_21,semag_conc_21]=conc_time_courses_mod(semag,2,tau,0,dos_switch_time);
[t_tirz5_21,tirz5_conc_21]=conc_time_courses_mod(tirz,2,tau,5,dos_switch_time);
[t_tirz10_21,tirz10_conc_21]=conc_time_courses_mod(tirz,2,tau,10,dos_switch_time);
[t_tirz15_21,tirz15_conc_21]=conc_time_courses_mod(tirz,2,tau,15,dos_switch_time);

tau=28*24;
[t_semag_28,semag_conc_28]=conc_time_courses_mod(semag,2,tau,0,dos_switch_time);
[t_tirz5_28,tirz5_conc_28]=conc_time_courses_mod(tirz,2,tau,5,dos_switch_time);
[t_tirz10_28,tirz10_conc_28]=conc_time_courses_mod(tirz,2,tau,10,dos_switch_time);
[t_tirz15_28,tirz15_conc_28]=conc_time_courses_mod(tirz,2,tau,15,dos_switch_time);


%%Placebo 1 
[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag_pl1,semag_conc_pl1,semag),t_semag_pl1,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(semag.bw0,y(:,3), bw_treat(semag,semag_conc_pl1), y(:,1),y(:,2));
semag_bw_pl1=100*(bw_treated-semag.bw0)/semag.bw0;


[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz5_pl1,tirz5_conc_pl1,tirz),t_tirz5_pl1,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz5_conc_pl1), y(:,1),y(:,2));
tirz5_bw_pl1=100*(bw_treated-tirz.bw0)/tirz.bw0;

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz10_pl1,tirz10_conc_pl1,tirz),t_tirz10_pl1,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz10_conc_pl1), y(:,1),y(:,2));
tirz10_bw_pl1=100*(bw_treated-tirz.bw0)/tirz.bw0;

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz15_pl1,tirz15_conc_pl1,tirz),t_tirz15_pl1,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz15_conc_pl1), y(:,1),y(:,2));
tirz15_bw_pl1=100*(bw_treated-tirz.bw0)/tirz.bw0;

%%Placebo 2 
[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag_pl2,semag_conc_pl2,semag),t_semag_pl2,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(semag.bw0,y(:,3), bw_treat(semag,semag_conc_pl2), y(:,1),y(:,2));
semag_bw_pl2=100*(bw_treated-semag.bw0)/semag.bw0;


[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz5_pl2,tirz5_conc_pl2,tirz),t_tirz5_pl2,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz5_conc_pl2), y(:,1),y(:,2));
tirz5_bw_pl2=100*(bw_treated-tirz.bw0)/tirz.bw0;

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz10_pl2,tirz10_conc_pl2,tirz),t_tirz10_pl2,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz10_conc_pl2), y(:,1),y(:,2));
tirz10_bw_pl2=100*(bw_treated-tirz.bw0)/tirz.bw0;

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz15_pl2,tirz15_conc_pl2,tirz),t_tirz15_pl2,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz15_conc_pl2), y(:,1),y(:,2));
tirz15_bw_pl2=100*(bw_treated-tirz.bw0)/tirz.bw0;


%%7-Day Dosing
[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag_7,semag_conc_7,semag),t_semag_7,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(semag.bw0,y(:,3), bw_treat(semag,semag_conc_7), y(:,1),y(:,2));
bw_treatedNOTRANSIENT=body_weightNOTRANSIENT(semag.bw0,y(:,3), bw_treat(semag,semag_conc_7), y(:,1),y(:,2));
semag_bw_7=100*(bw_treated-semag.bw0)/semag.bw0;
semag_bw_7NOTRANSIENT=100*(bw_treatedNOTRANSIENT-semag.bw0)/semag.bw0;

semag_bw_7_av=filter(ones(1,7*24)/(7*24),1,semag_bw_7);
semag_7_av_ss=semag_bw_7_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz5_7,tirz5_conc_7,tirz),t_tirz5_7,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz5_conc_7), y(:,1),y(:,2));
tirz5_bw_7=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz5_bw_7_av=filter(ones(1,7*24)/(7*24),1,tirz5_bw_7);
tirz5_7_av_ss=tirz5_bw_7_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz10_7,tirz10_conc_7,tirz),t_tirz10_7,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz10_conc_7), y(:,1),y(:,2));
tirz10_bw_7=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz10_bw_7_av=filter(ones(1,7*24)/(7*24),1,tirz10_bw_7);
tirz10_7_av_ss=tirz10_bw_7_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz15_7,tirz15_conc_7,tirz),t_tirz15_7,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz15_conc_7), y(:,1),y(:,2));
tirz15_bw_7=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz15_bw_7_av=filter(ones(1,7*24)/(7*24),1,tirz15_bw_7);
tirz15_7_av_ss=tirz15_bw_7_av(tpoint);

%%10-Day Dosing
[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag_10,semag_conc_10,semag),t_semag_10,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(semag.bw0,y(:,3), bw_treat(semag,semag_conc_10), y(:,1),y(:,2));
semag_bw_10=100*(bw_treated-semag.bw0)/semag.bw0;

semag_bw_10_av=filter(ones(1,7*24)/(7*24),1,semag_bw_10);
semag_10_av_ss=semag_bw_10_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz5_10,tirz5_conc_10,tirz),t_tirz5_10,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz5_conc_10), y(:,1),y(:,2));
tirz5_bw_10=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz5_bw_10_av=filter(ones(1,7*24)/(7*24),1,tirz5_bw_10);
tirz5_10_av_ss=tirz5_bw_10_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz10_10,tirz10_conc_10,tirz),t_tirz10_10,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz10_conc_10), y(:,1),y(:,2));
tirz10_bw_10=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz10_bw_10_av=filter(ones(1,7*24)/(7*24),1,tirz10_bw_10);
tirz10_10_av_ss=tirz10_bw_10_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz15_10,tirz15_conc_10,tirz),t_tirz15_10,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz15_conc_10), y(:,1),y(:,2));
tirz15_bw_10=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz15_bw_10_av=filter(ones(1,7*24)/(7*24),1,tirz15_bw_10);
tirz15_10_av_ss=tirz15_bw_10_av(tpoint);

%%14-Day Dosing
[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag_14,semag_conc_14,semag),t_semag_14,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(semag.bw0,y(:,3), bw_treat(semag,semag_conc_14), y(:,1),y(:,2));
semag_bw_14=100*(bw_treated-semag.bw0)/semag.bw0;

semag_bw_14_av=filter(ones(1,7*24)/(7*24),1,semag_bw_14);
semag_14_av_ss=semag_bw_14_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz5_14,tirz5_conc_14,tirz),t_tirz5_14,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz5_conc_14), y(:,1),y(:,2));
tirz5_bw_14=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz5_bw_14_av=filter(ones(1,7*24)/(7*24),1,tirz5_bw_14);
tirz5_14_av_ss=tirz5_bw_14_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz10_14,tirz10_conc_14,tirz),t_tirz10_14,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz10_conc_14), y(:,1),y(:,2));
tirz10_bw_14=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz10_bw_14_av=filter(ones(1,7*24)/(7*24),1,tirz10_bw_14);
tirz10_14_av_ss=tirz10_bw_14_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz15_14,tirz15_conc_14,tirz),t_tirz15_14,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz15_conc_14), y(:,1),y(:,2));
tirz15_bw_14=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz15_bw_14_av=filter(ones(1,7*24)/(7*24),1,tirz15_bw_14);
tirz15_14_av_ss=tirz15_bw_14_av(tpoint);

%%21-Day Dosing
[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag_21,semag_conc_21,semag),t_semag_21,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(semag.bw0,y(:,3), bw_treat(semag,semag_conc_21), y(:,1),y(:,2));
semag_bw_21=100*(bw_treated-semag.bw0)/semag.bw0;

semag_bw_21_av=filter(ones(1,7*24)/(7*24),1,semag_bw_21);
semag_21_av_ss=semag_bw_21_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz5_21,tirz5_conc_21,tirz),t_tirz5_21,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz5_conc_21), y(:,1),y(:,2));
tirz5_bw_21=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz5_bw_21_av=filter(ones(1,7*24)/(7*24),1,tirz5_bw_21);
tirz5_21_av_ss=tirz5_bw_21_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz10_21,tirz10_conc_21,tirz),t_tirz10_21,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz10_conc_21), y(:,1),y(:,2));
tirz10_bw_21=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz10_bw_21_av=filter(ones(1,7*24)/(7*24),1,tirz10_bw_21);
tirz10_21_av_ss=tirz10_bw_21_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz15_21,tirz15_conc_21,tirz),t_tirz15_21,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz15_conc_21), y(:,1),y(:,2));
tirz15_bw_21=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz15_bw_21_av=filter(ones(1,7*24)/(7*24),1,tirz15_bw_21);
tirz15_21_av_ss=tirz15_bw_21_av(tpoint);

%%28-Day Dosing
[~,y]=ode45(@(t,y) ode_solver(t,y,t_semag_28,semag_conc_28,semag),t_semag_28,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(semag.bw0,y(:,3), bw_treat(semag,semag_conc_28), y(:,1),y(:,2));
semag_bw_28=100*(bw_treated-semag.bw0)/semag.bw0;

semag_bw_28_av=filter(ones(1,7*24)/(7*24),1,semag_bw_28);
semag_28_av_ss=semag_bw_28_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz5_28,tirz5_conc_28,tirz),t_tirz5_28,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz5_conc_28), y(:,1),y(:,2));
tirz5_bw_28=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz5_bw_28_av=filter(ones(1,7*24)/(7*24),1,tirz5_bw_28);
tirz5_28_av_ss=tirz5_bw_28_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz10_28,tirz10_conc_28,tirz),t_tirz10_28,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz10_conc_28), y(:,1),y(:,2));
tirz10_bw_28=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz10_bw_28_av=filter(ones(1,7*24)/(7*24),1,tirz10_bw_28);
tirz10_28_av_ss=tirz10_bw_28_av(tpoint);

[~,y]=ode45(@(t,y) ode_solver(t,y,t_tirz15_28,tirz15_conc_28,tirz),t_tirz15_28,[bw_pld_ic;bw_pln_ic;bw_s_ic],ode_options);
bw_treated=body_weight(tirz.bw0,y(:,3), bw_treat(tirz,tirz15_conc_28), y(:,1),y(:,2));
tirz15_bw_28=100*(bw_treated-tirz.bw0)/tirz.bw0;

tirz15_bw_28_av=filter(ones(1,7*24)/(7*24),1,tirz15_bw_28);
tirz15_28_av_ss=tirz15_bw_28_av(tpoint);

%% For table
av_ss=[semag_bw_pl1(tpoint), semag_bw_pl2(tpoint),semag_7_av_ss,semag_10_av_ss,semag_14_av_ss,semag_21_av_ss,semag_28_av_ss; ...
tirz5_bw_pl1(tpoint), tirz5_bw_pl2(tpoint),tirz5_7_av_ss,tirz5_10_av_ss,tirz5_14_av_ss,tirz5_21_av_ss,tirz5_28_av_ss; ...
tirz10_bw_pl1(tpoint), tirz10_bw_pl2(tpoint),tirz10_7_av_ss,tirz10_10_av_ss,tirz10_14_av_ss,tirz10_21_av_ss,tirz10_28_av_ss;...
tirz15_bw_pl1(tpoint), tirz15_bw_pl2(tpoint),tirz15_7_av_ss,tirz15_10_av_ss,tirz15_14_av_ss,tirz15_21_av_ss,tirz15_28_av_ss];

av_ss=round(av_ss)

%% save data

% timestamp = char(datetime('now'));
% timestamp = strrep(timestamp,' ','-');
% timestamp = strrep(timestamp,':','');
% filename = strcat(['data_Case_',timestamp,'.mat']);
filename = strcat(['data_Case.mat']);
save(filename);

%% Graphing

col_vec=["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628"];
mar_vec=["diamond","v","square","o","hexagram"];

linwid=1.5;
m=4;

hfig=figure()
subplot(2,2,1)

plot(t_semag_pl1/(7*24),semag_bw_pl1,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','--')
hold on;
plot(t_semag_pl2/(7*24),semag_bw_pl2,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','-.')
plot(t_semag_7/(7*24),semag_bw_7,'LineWidth', linwid,'Color',col_vec(2))
plot(t_semag_10/(7*24),semag_bw_10,'LineWidth',linwid,'Color',col_vec(3))
plot(t_semag_14/(7*24),semag_bw_14,'LineWidth',linwid,'Color',col_vec(4))
plot(t_semag_21/(7*24),semag_bw_21,'LineWidth',linwid,'Color',col_vec(5))
plot(t_semag_28/(7*24),semag_bw_28,'LineWidth',linwid,'Color',col_vec(6))
ylabel("Percent Change in BW")
ylim([-25,2])
yline(0,'LineStyle',":",'LineWidth',linwid)
grid on 
grid minor
box on 

xlabel("Time (weeks)")
xlim([0,tpoint/(7*24)])
xticks([0,40,80,120,160,200,240])
text(65,0.85,"Baseline BW=104.8 kg")
yline(0,'LineStyle',":",'LineWidth',linwid)
xlt=title("Semaglutide (2.4 mg) Weight Loss")

grid on 
grid minor


subplot(2,2,2)
plot(t_tirz5_pl1/(7*24),tirz5_bw_pl1,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','--')
hold on
plot(t_tirz5_pl2/(7*24),tirz5_bw_pl2,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','-.')
plot(t_tirz5_7/(7*24),tirz5_bw_7,'LineWidth', linwid,'Color',col_vec(2))
plot(t_tirz5_10/(7*24),tirz5_bw_10,'LineWidth',linwid,'Color',col_vec(3))
plot(t_tirz5_14/(7*24),tirz5_bw_14,'LineWidth',linwid,'Color',col_vec(4))
plot(t_tirz5_21/(7*24),tirz5_bw_21,'LineWidth',linwid,'Color',col_vec(5))
plot(t_tirz5_28/(7*24),tirz5_bw_28,'LineWidth',linwid,'Color',col_vec(6))
yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
ylabel("Percent Change in BW")

xticks([0,40,80,120,160,200,240])
xlim([0,tpoint/(7*24)])
ylim([-25,2])

xlt=title("Tirzepatide (5 mg)  Weight Loss")
text(65,0.85,"Baseline BW=104.8 kg")
yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
xlh=xlabel("Time (weeks)")

subplot(2,2,3)
plot(t_tirz10_pl1/(7*24),tirz10_bw_pl1,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','--')

hold on
plot(t_tirz10_pl2/(7*24),tirz10_bw_pl2,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','-.')
plot(t_tirz10_7/(7*24),tirz10_bw_7,'LineWidth', linwid,'Color',col_vec(2))
plot(t_tirz10_10/(7*24),tirz10_bw_10,'LineWidth',linwid,'Color',col_vec(3))
plot(t_tirz10_14/(7*24),tirz10_bw_14,'LineWidth',linwid,'Color',col_vec(4))
plot(t_tirz10_21/(7*24),tirz10_bw_21,'LineWidth',linwid,'Color',col_vec(5))
plot(t_tirz10_28/(7*24),tirz10_bw_28,'LineWidth',linwid,'Color',col_vec(6))
yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
ylabel("Percent Change in BW")

xticks([0,40,80,120,160,200,240])
xlim([0,tpoint/(7*24)])
ylim([-25,2])
xlt=title("Tirzepatide (10 mg)  Weight Loss")
text(65,0.85,"Baseline BW=104.8 kg")

yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
xlh=xlabel("Time (weeks)")

subplot(2,2,4)
plot(t_tirz15_pl1/(7*24),tirz15_bw_pl1,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','--')
hold on
plot(t_tirz15_pl2/(7*24),tirz15_bw_pl2,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','-.')
plot(t_tirz15_7/(7*24),tirz15_bw_7,'LineWidth', linwid,'Color',col_vec(2))
plot(t_tirz15_10/(7*24),tirz15_bw_10,'LineWidth',linwid,'Color',col_vec(3))
plot(t_tirz15_14/(7*24),tirz15_bw_14,'LineWidth',linwid,'Color',col_vec(4))
plot(t_tirz15_21/(7*24),tirz15_bw_21,'LineWidth',linwid,'Color',col_vec(5))
plot(t_tirz15_28/(7*24),tirz15_bw_28,'LineWidth',linwid,'Color',col_vec(6))
yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
ylabel("Percent Change in BW")

xticks([0,40,80,120,160,200,240])
xlim([0,tpoint/(7*24)])
ylim([-25,2])

xlt=title("Tirzepatide (15 mg)  Weight Loss")
text(65,0.85,"Baseline BW=104.8 kg")
yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
xlh=xlabel("Time (weeks)")

lgd=legend("$\textup{PBO}_{1}$","$\textup{PBO}_{2}$","$\tau=7\textup{ days}$","$\tau=10\textup{ days}$", "$\tau=14\textup{ days}$","$\tau=21\textup{ days}$","$\tau=28\textup{ days}$",'NumColumns',7,'Location','south')

fname="casereport_graph";
picturewidth = 40; % set this parameter and keep it forever
hw_ratio = 1; % feel free to play with this ratio
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

%%

function bw=body_weight(bw_0, bw_s, bw_i, bw_pld, bw_pln)
bw=bw_0-bw_s-bw_i+bw_pld+bw_pln;
end

function bw=body_weightNOTRANSIENT(bw_0, bw_s, bw_i, bw_pld, bw_pln)
bw=bw_0-bw_s-bw_i+bw_pln;
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

function bwi=bw_treat(par,cavg)
%bwi=0.01*bw_0*emaxi*cavg/(ec50+cavg);
bwi=par.bw0*par.emaxi*cavg./(par.ec50i+cavg);
end



