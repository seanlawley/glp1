close all
clear all

load('data_Case.mat')
close all

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


%% Graphing

col_vec=["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628"];
mar_vec=["diamond","v","square","o","hexagram"];

linwid=1.5;
m=4;

hfig=figure()
tcl = tiledlayout(2,2,'TileSpacing', 'loose');

%%%%%
% subplot(2,2,1)
nexttile(tcl)

plot(t_semag_pl1/(7*24),semag_bw_pl1,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','--')
hold on;
plot(t_semag_pl2/(7*24),semag_bw_pl2,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','-.')
plot(t_semag_10/(7*24),semag_bw_10,'LineWidth',linwid,'Color',col_vec(3))
plot(t_semag_14/(7*24),semag_bw_14,'LineWidth',linwid,'Color',col_vec(4))
plot(t_semag_21/(7*24),semag_bw_21,'LineWidth',linwid,'Color',col_vec(5))
plot(t_semag_28/(7*24),semag_bw_28,'LineWidth',linwid,'Color',col_vec(6))
plot(t_semag_7/(7*24),semag_bw_7,'LineWidth', linwid,'Color',col_vec(2))

e1=errorbar(t_pl_data,pl_data,pl_data-pl_bot,pl_top-pl_data,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');
e1.Color='black';
e1.LineWidth=linwid;
e2=errorbar(t_semag_data,semag_data,semag_data-semag_bot,semag_top-semag_data,'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');
e2.Color='black';
e2.LineWidth=linwid;

semag_7mid=trapz(t_semag_7(7*24*120+1:7*24*121+1),semag_bw_7(7*24*120+1:7*24*121+1))/(7*24);
semag_pl1ss=trapz(t_semag_pl1(7*24*240+1:7*24*241+1),semag_bw_pl1(7*24*240+1:7*24*241+1))/(7*24);
semag_7ss=trapz(t_semag_7(7*24*240+1:7*24*241+1),semag_bw_7(7*24*240+1:7*24*241+1))/(7*24);
semag_10ss=trapz(t_semag_10(7*24*240+1:7*24*241+1),semag_bw_10(7*24*240+1:7*24*241+1))/(7*24);
semag_14ss=trapz(t_semag_14(7*24*240+1:7*24*241+1),semag_bw_14(7*24*240+1:7*24*241+1))/(7*24);
semag_21ss=trapz(t_semag_21(7*24*240+1:7*24*241+1),semag_bw_21(7*24*240+1:7*24*241+1))/(7*24);
semag_28ss=trapz(t_semag_28(7*24*240+1:7*24*241+1),semag_bw_28(7*24*240+1:7*24*241+1))/(7*24);
text(114, semag_7mid-1,sprintf('%.0f', semag_7mid), 'FontSize', 35, 'Color', col_vec(2), 'HorizontalAlignment', 'left');
text(241, semag_pl1ss, sprintf('%.0f', semag_pl1ss), 'FontSize', 35, 'Color', col_vec(1), 'HorizontalAlignment', 'left');
text(241, semag_7ss, sprintf('%.0f', semag_7ss), 'FontSize', 35, 'Color', col_vec(2), 'HorizontalAlignment', 'left');
text(241, semag_10ss, sprintf('%.0f', semag_10ss), 'FontSize', 35, 'Color', col_vec(3), 'HorizontalAlignment', 'left');
text(241, semag_14ss, sprintf('%.0f', semag_14ss), 'FontSize', 35, 'Color', col_vec(4), 'HorizontalAlignment', 'left');
text(241, semag_21ss, sprintf('%.0f', semag_21ss), 'FontSize', 35, 'Color', col_vec(5), 'HorizontalAlignment', 'left');
text(241, semag_28ss, sprintf('%.0f', semag_28ss), 'FontSize', 35, 'Color', col_vec(6), 'HorizontalAlignment', 'left');


ylabel("Percent Change in BW")
ylim([-26,0])
grid on 
grid minor
box on 

xlabel("Time (weeks)")
xlim([0,tpoint/(7*24)])
xticks([0,40,80,120,160,200,240])
% text(65,0.85,"Baseline BW=104.8 kg")
% yline(0,'LineStyle',":",'LineWidth',linwid)
xlt=title("Semaglutide (2.4 mg) Weight Loss")

grid on 
grid minor



%%%%%
load('dataJM.mat')
% subplot(2,2,2)
nexttile(tcl)

plot(t_tirz5_pl1/(7*24),tirz5_bw_pl1,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','--')
hold on
plot(t_tirz5_pl2/(7*24),tirz5_bw_pl2,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','-.')
plot(t_tirz5_10/(7*24),tirz5_bw_10,'LineWidth',linwid,'Color',col_vec(3))
plot(t_tirz5_14/(7*24),tirz5_bw_14,'LineWidth',linwid,'Color',col_vec(4))
plot(t_tirz5_21/(7*24),tirz5_bw_21,'LineWidth',linwid,'Color',col_vec(5))
plot(t_tirz5_28/(7*24),tirz5_bw_28,'LineWidth',linwid,'Color',col_vec(6))
plot(t_tirz5_7/(7*24),tirz5_bw_7,'LineWidth', linwid,'Color',col_vec(2))

e1=errorbar(JM(:,1),JM(:,2),JM(:,2)-pl_bot,pl_top-JM(:,2),'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');
e1.Color='black';
e1.LineWidth=linwid;
e2=errorbar(JM(:,7),JM(:,8),JM(:,8)-t_5_bot,t_5_top-JM(:,8),'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');
e2.Color='black';
e2.LineWidth=linwid;

grid on 
grid minor

tirz5_7mid=trapz(t_tirz5_7(7*24*120+1:7*24*121+1),tirz5_bw_7(7*24*120+1:7*24*121+1))/(7*24);
tirz5_pl1ss=trapz(t_tirz5_pl1(7*24*240+1:7*24*241+1),tirz5_bw_pl1(7*24*240+1:7*24*241+1))/(7*24);
tirz5_pl2ss=trapz(t_tirz5_pl2(7*24*240+1:7*24*241+1),tirz5_bw_pl2(7*24*240+1:7*24*241+1))/(7*24);
tirz5_7ss=trapz(t_tirz5_7(7*24*240+1:7*24*241+1),tirz5_bw_7(7*24*240+1:7*24*241+1))/(7*24);
tirz5_10ss=trapz(t_tirz5_10(7*24*240+1:7*24*241+1),tirz5_bw_10(7*24*240+1:7*24*241+1))/(7*24);
tirz5_14ss=trapz(t_tirz5_14(7*24*240+1:7*24*241+1),tirz5_bw_14(7*24*240+1:7*24*241+1))/(7*24);
tirz5_21ss=trapz(t_tirz5_21(7*24*240+1:7*24*241+1),tirz5_bw_21(7*24*240+1:7*24*241+1))/(7*24);
tirz5_28ss=trapz(t_tirz5_28(7*24*240+1:7*24*241+1),tirz5_bw_28(7*24*240+1:7*24*241+1))/(7*24);
text(114, tirz5_7mid-1,sprintf('%.0f', tirz5_7mid), 'FontSize', 35, 'Color', col_vec(2), 'HorizontalAlignment', 'left');
text(241, tirz5_pl1ss, sprintf('%.0f', tirz5_pl1ss), 'FontSize', 35, 'Color', col_vec(1), 'HorizontalAlignment', 'left');
text(241, tirz5_7ss, sprintf('%.0f', tirz5_7ss), 'FontSize', 35, 'Color', col_vec(2), 'HorizontalAlignment', 'left');
text(241, tirz5_10ss, sprintf('%.0f', tirz5_10ss), 'FontSize', 35, 'Color', col_vec(3), 'HorizontalAlignment', 'left');
text(241, tirz5_14ss, sprintf('%.0f', tirz5_14ss), 'FontSize', 35, 'Color', col_vec(4), 'HorizontalAlignment', 'left');
text(241, tirz5_21ss, sprintf('%.0f', tirz5_21ss), 'FontSize', 35, 'Color', col_vec(5), 'HorizontalAlignment', 'left');
text(241, tirz5_28ss, sprintf('%.0f', tirz5_28ss), 'FontSize', 35, 'Color', col_vec(6), 'HorizontalAlignment', 'left');


ylabel("Percent Change in BW")

xticks([0,40,80,120,160,200,240])
xlim([0,tpoint/(7*24)])
ylim([-26,0])

xlt=title("Tirzepatide (5 mg)  Weight Loss")
% text(65,0.85,"Baseline BW=104.8 kg")
% yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
xlh=xlabel("Time (weeks)")


%%%%%
% subplot(2,2,3)
nexttile(tcl)

plot(t_tirz10_pl1/(7*24),tirz10_bw_pl1,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','--')

hold on
plot(t_tirz10_pl2/(7*24),tirz10_bw_pl2,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','-.')
plot(t_tirz10_10/(7*24),tirz10_bw_10,'LineWidth',linwid,'Color',col_vec(3))
plot(t_tirz10_14/(7*24),tirz10_bw_14,'LineWidth',linwid,'Color',col_vec(4))
plot(t_tirz10_21/(7*24),tirz10_bw_21,'LineWidth',linwid,'Color',col_vec(5))
plot(t_tirz10_28/(7*24),tirz10_bw_28,'LineWidth',linwid,'Color',col_vec(6))
plot(t_tirz10_7/(7*24),tirz10_bw_7,'LineWidth', linwid,'Color',col_vec(2))

e1=errorbar(JM(:,1),JM(:,2),JM(:,2)-pl_bot,pl_top-JM(:,2),'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');
e1.Color='black';
e1.LineWidth=linwid;
e3=errorbar(JM(:,13),JM(:,14),JM(:,14)-t_10_bot,t_10_top-JM(:,14),'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');
e3.Color='black';
e3.LineWidth=linwid;


grid on 
grid minor

tirz10_7mid=trapz(t_tirz10_7(7*24*120+1:7*24*121+1),tirz10_bw_7(7*24*120+1:7*24*121+1))/(7*24);
tirz10_pl1ss=trapz(t_tirz10_pl1(7*24*240+1:7*24*241+1),tirz10_bw_pl1(7*24*240+1:7*24*241+1))/(7*24);
tirz10_7ss=trapz(t_tirz10_7(7*24*240+1:7*24*241+1),tirz10_bw_7(7*24*240+1:7*24*241+1))/(7*24);
tirz10_10ss=trapz(t_tirz10_10(7*24*240+1:7*24*241+1),tirz10_bw_10(7*24*240+1:7*24*241+1))/(7*24);
tirz10_14ss=trapz(t_tirz10_14(7*24*240+1:7*24*241+1),tirz10_bw_14(7*24*240+1:7*24*241+1))/(7*24);
tirz10_21ss=trapz(t_tirz10_21(7*24*240+1:7*24*241+1),tirz10_bw_21(7*24*240+1:7*24*241+1))/(7*24);
tirz10_28ss=trapz(t_tirz10_28(7*24*240+1:7*24*241+1),tirz10_bw_28(7*24*240+1:7*24*241+1))/(7*24);
text(114, tirz10_7mid-1,sprintf('%.0f', tirz10_7mid), 'FontSize', 35, 'Color', col_vec(2), 'HorizontalAlignment', 'left');
text(241, tirz10_pl1ss, sprintf('%.0f', tirz10_pl1ss), 'FontSize', 35, 'Color', col_vec(1), 'HorizontalAlignment', 'left');
text(241, tirz10_7ss, sprintf('%.0f', tirz10_7ss), 'FontSize', 35, 'Color', col_vec(2), 'HorizontalAlignment', 'left');
text(241, tirz10_10ss, sprintf('%.0f', tirz10_10ss), 'FontSize', 35, 'Color', col_vec(3), 'HorizontalAlignment', 'left');
text(241, tirz10_14ss, sprintf('%.0f', tirz10_14ss), 'FontSize', 35, 'Color', col_vec(4), 'HorizontalAlignment', 'left');
text(241, tirz10_21ss, sprintf('%.0f', tirz10_21ss), 'FontSize', 35, 'Color', col_vec(5), 'HorizontalAlignment', 'left');
text(241, tirz10_28ss, sprintf('%.0f', tirz10_28ss), 'FontSize', 35, 'Color', col_vec(6), 'HorizontalAlignment', 'left');


ylabel("Percent Change in BW")

xticks([0,40,80,120,160,200,240])
xlim([0,tpoint/(7*24)])
ylim([-26,0])
xlt=title("Tirzepatide (10 mg)  Weight Loss")
% text(65,0.85,"Baseline BW=104.8 kg")
% yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
xlh=xlabel("Time (weeks)")

%%%%%
% subplot(2,2,4)
nexttile(tcl)

p1=plot(t_tirz15_pl1/(7*24),tirz15_bw_pl1,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','--');
hold on
p2=plot(t_tirz15_pl2/(7*24),tirz15_bw_pl2,'LineWidth',linwid,'Color',col_vec(1),'LineStyle','-.');
p4=plot(t_tirz15_10/(7*24),tirz15_bw_10,'LineWidth',linwid,'Color',col_vec(3));
p5=plot(t_tirz15_14/(7*24),tirz15_bw_14,'LineWidth',linwid,'Color',col_vec(4));
p6=plot(t_tirz15_21/(7*24),tirz15_bw_21,'LineWidth',linwid,'Color',col_vec(5));
p7=plot(t_tirz15_28/(7*24),tirz15_bw_28,'LineWidth',linwid,'Color',col_vec(6));
p3=plot(t_tirz15_7/(7*24),tirz15_bw_7,'LineWidth', linwid,'Color',col_vec(2));

e1=errorbar(JM(:,1),JM(:,2),JM(:,2)-pl_bot,pl_top-JM(:,2),'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');
e1.Color='black';
e1.LineWidth=linwid;
e4=errorbar(JM(:,19),JM(:,20),JM(:,20)-t_15_bot,t_15_top-JM(:,20),'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');
e4.Color='black';
e4.LineWidth=linwid;


grid on 
grid minor

tirz15_7mid=trapz(t_tirz15_7(7*24*120+1:7*24*121+1),tirz15_bw_7(7*24*120+1:7*24*121+1))/(7*24);
tirz15_pl1ss=trapz(t_tirz15_pl1(7*24*240+1:7*24*241+1),tirz15_bw_pl1(7*24*240+1:7*24*241+1))/(7*24);
tirz15_7ss=trapz(t_tirz15_7(7*24*240+1:7*24*241+1),tirz15_bw_7(7*24*240+1:7*24*241+1))/(7*24);
tirz15_10ss=trapz(t_tirz15_10(7*24*240+1:7*24*241+1),tirz15_bw_10(7*24*240+1:7*24*241+1))/(7*24);
tirz15_14ss=trapz(t_tirz15_14(7*24*240+1:7*24*241+1),tirz15_bw_14(7*24*240+1:7*24*241+1))/(7*24);
tirz15_21ss=trapz(t_tirz15_21(7*24*240+1:7*24*241+1),tirz15_bw_21(7*24*240+1:7*24*241+1))/(7*24);
tirz15_28ss=trapz(t_tirz15_28(7*24*240+1:7*24*241+1),tirz15_bw_28(7*24*240+1:7*24*241+1))/(7*24);
text(114, tirz15_7mid-1,sprintf('%.0f', tirz15_7mid), 'FontSize', 35, 'Color', col_vec(2), 'HorizontalAlignment', 'left');
text(241, tirz15_pl1ss, sprintf('%.0f', tirz15_pl1ss), 'FontSize', 35, 'Color', col_vec(1), 'HorizontalAlignment', 'left');
text(241, tirz15_7ss, sprintf('%.0f', tirz15_7ss), 'FontSize', 35, 'Color', col_vec(2), 'HorizontalAlignment', 'left');
text(241, tirz15_10ss, sprintf('%.0f', tirz15_10ss), 'FontSize', 35, 'Color', col_vec(3), 'HorizontalAlignment', 'left');
text(241, tirz15_14ss, sprintf('%.0f', tirz15_14ss), 'FontSize', 35, 'Color', col_vec(4), 'HorizontalAlignment', 'left');
text(241, tirz15_21ss, sprintf('%.0f', tirz15_21ss), 'FontSize', 35, 'Color', col_vec(5), 'HorizontalAlignment', 'left');
text(241, tirz15_28ss, sprintf('%.0f', tirz15_28ss), 'FontSize', 35, 'Color', col_vec(6), 'HorizontalAlignment', 'left');


ylabel("Percent Change in BW")

xticks([0,40,80,120,160,200,240])
xlim([0,tpoint/(7*24)])
ylim([-26,0])

xlt=title("Tirzepatide (15 mg)  Weight Loss")
% text(65,0.85,"Baseline BW=104.8 kg")
% yline(0,'LineStyle',":",'LineWidth',linwid)%legend("Placebo","Placebo (Post 24 Weeks)","15mg, DI of 7 Days","15mg, DI of 10 Days", "15mg, DI of 14 Days","15mg, DI of 21 Days",'Location','southwest')
grid on 
grid minor
xlh=xlabel("Time (weeks)")

lgd=legend([p1,p2,p3,p4,p5,p6,p7],"$\textup{PBO}_{1}$","$\textup{PBO}_{2}$","$\tau=7\textup{ days}$","$\tau=10\textup{ days}$", "$\tau=14\textup{ days}$","$\tau=21\textup{ days}$","$\tau=28\textup{ days}$",'NumColumns',7);

% Construct a Legend with the data from the sub-plots
% hL = legend([line1,line2,line3,line4]); 
% Move the legend to the right side of the figure
lgd.Layout.Tile = 'South';

%% display table

table4word=[7/10 7/14 7/21 7/28;...
    semag_10ss/semag_7ss semag_14ss/semag_7ss semag_21ss/semag_7ss semag_28ss/semag_7ss;...
    tirz5_10ss/tirz5_7ss tirz5_14ss/tirz5_7ss tirz5_21ss/tirz5_7ss tirz5_28ss/tirz5_7ss;...
    tirz10_10ss/tirz10_7ss tirz10_14ss/tirz10_7ss tirz10_21ss/tirz10_7ss tirz10_28ss/tirz10_7ss;...
    tirz15_10ss/tirz15_7ss tirz15_14ss/tirz15_7ss tirz15_21ss/tirz15_7ss tirz15_28ss/tirz15_7ss;...
    ];

table4word=round(100*table4word)


%% save and format figure

fname="figCase";
picturewidth = 40; % set this parameter and keep it forever
hw_ratio = .8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')

xDist=xlim;
xDist=xDist(end)-xDist(1);
yDist=ylim;
yDist=yDist(end)-yDist(1);
text(-1.375*xDist,1.36*yDist, 'A', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(-.13*xDist,1.36*yDist, 'B', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(-1.375*xDist,.05*yDist, 'C', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(-.13*xDist,.05*yDist, 'D', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');

legend boxoff
%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(hfig,strcat(fname,'.fig'))