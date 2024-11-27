%% load data
close all
clear all
load('data_SteadyState')
close all
%% plotting

linwid=1.5;
hfig=figure()
hold on;

%% panel A
subplot(1,2,1)
hold all

plot(dosing_interval/(24),semag_av,'-','Color',col_vec(1),'LineWidth',2);
% plot(dosing_interval/(24),tirz_av_5,'--','Color',col_vec(2),'LineWidth',2);
% plot(dosing_interval/(24),tirz_av_10,':','Color',col_vec(3),'LineWidth',2);
% plot(dosing_interval/(24),tirz_av_15,'-.','Color',col_vec(4),'LineWidth',2);

ms=10;
lw=1.5;
plot([7 14],[semag_av(1) semag_av(8)],'s','Color',col_vec(1),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(1))
% plot([7 14],[tirz_av_5(1) tirz_av_5(8)],'o','Color',col_vec(2),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(2))
% plot([7 14],[tirz_av_10(1) tirz_av_10(8)],'v','Color',col_vec(3),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(3))
% plot([7 14],[tirz_av_15(1) tirz_av_15(8)],'diamond','Color',col_vec(4),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(4))
legend("Semaglutide 2.4 mg",...
    "Location","Southeast")

xlabel("Dosing Interval (Days)")
ylabel("Percent Change in Body Weight")
grid on 
grid minor
% xlim([4,29])
xlim([6.5,28.5])
ylim([-20,0])
%ylabel("$\%\Delta$ in BW ($\%$)")
legend boxoff


xticks([7,10,14,17,21,24,28])
% xticklabels(['7','14','10','17','21','24','28'])

h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.XAxis.MinorTickValues = 7:1:28; % Minor ticks which don't line up with majors

%% panel B
subplot(1,2,2)
hold all


% scatter(dosing_interval/(24),semag_av,80,'MarkerEdgeColor',col_vec(1),'MarkerFaceColor',col_vec(1),'Marker',"square");
% scatter(dosing_interval/(24),tirz_av_5,80,'MarkerEdgeColor',col_vec(2),'MarkerFaceColor',col_vec(2),'Marker',"o");
% scatter(dosing_interval/(24),tirz_av_10,80,'MarkerEdgeColor',col_vec(3),'MarkerFaceColor',col_vec(3),'Marker',"v");
% scatter(dosing_interval/(24),tirz_av_15,80,'MarkerEdgeColor',col_vec(4),'MarkerFaceColor',col_vec(4),'Marker',"^");

cost=dosing_interval/(24);
cost=cost(1)./cost;

plot(cost,semag_av./semag_av(1),'-','Color',col_vec(1),'LineWidth',2);
% plot(cost,tirz_av_5./tirz_av_5(1),'--','Color',col_vec(2),'LineWidth',2);
% plot(cost,tirz_av_10./tirz_av_10(1),':','Color',col_vec(3),'LineWidth',2);
% plot(cost,tirz_av_15./tirz_av_15(1),'-.','Color',col_vec(4),'LineWidth',2);
plot(cost,cost,'k','LineWidth',1)
% plot(dosing_interval/(24),semag_av./semag_av(1),'-','Color',col_vec(1),'LineWidth',2);
% plot(dosing_interval/(24),tirz_av_5./tirz_av_5(1),'--','Color',col_vec(2),'LineWidth',2);
% plot(dosing_interval/(24),tirz_av_10./tirz_av_10(1),':','Color',col_vec(3),'LineWidth',2);
% plot(dosing_interval/(24),tirz_av_15./tirz_av_15(1),'-.','Color',col_vec(4),'LineWidth',2);

%plot(dosing_interval/24,semag_ss*100,'LineWidth',linwid,'Color',col_vec(1));

%plot(dosing_interval/24,tirz_5*100,'LineWidth',linwid,'Color',col_vec(2));

%plot(dosing_interval/24,tirz_10*100,'LineWidth',linwid,'Color',col_vec(3));

%plot(dosing_interval/24,tirz_15*100,'LineWidth',linwid,'Color',col_vec(4));

ms=10;
lw=1.5;
plot([.5],[semag_av(8)]./semag_av(1),'s','Color',col_vec(1),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(1))
% plot([1 .5],[tirz_av_5(1) tirz_av_5(8)]./tirz_av_5(1),'o','Color',col_vec(2),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(2))
% plot([1 .5],[tirz_av_10(1) tirz_av_10(8)]./tirz_av_10(1),'v','Color',col_vec(3),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(3))
% plot([1 .5],[tirz_av_15(1) tirz_av_15(8)]./tirz_av_15(1),'diamond','Color',col_vec(4),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(4))
% plot([1 .5],[1 .5],'*k','MarkerSize',ms,'LineWidth',lw)

legend("Semaglutide 2.4 mg",...
    "Linear",...
    "Location","Southeast")

xlabel("Cost (relative to once-weekly)")
ylabel("Efficacy (relative to once-weekly)")
grid on 
grid minor
xlim([.25,1])
ylim([.25,1])
%ylabel("$\%\Delta$ in BW ($\%$)")
legend boxoff

text(.41, semag_av(8)./semag_av(1),  sprintf('%.0f\\%%', 100*semag_av(8)./semag_av(1)), 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');
% text(.42, .76, '73\% efficacy for 50\% cost', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');

xticks([.25,.5,.7,1])
yticks([.25,.5,.7,1])

h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
temp=flip(7./[7:1:28]);
h.XAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors
h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.YAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors


%% format and save figure
fname="figSemag"

picturewidth = 40; % set this parameter and keep it forever
hw_ratio = 0.4; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document


% set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')


text(-.85, 1, 'A', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(0.125, 1, 'B', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');

legend boxoff
%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(hfig,strcat(fname,'.fig'))