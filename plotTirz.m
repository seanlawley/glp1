%% load data
close all
clear all
load('data_SteadyState.mat')

%% Graphing

linwid=1.5;
hfig=figure()
hold on;

col_vec=["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33"];

%% panel A
subplot(2,2,1)
hold all

p1=plot(dosing_interval/(24),tirz_av_5,'--','Color',col_vec(2),'LineWidth',2);
p2=plot(dosing_interval/(24),tirz_av_10,':','Color',col_vec(3),'LineWidth',2);
p3=plot(dosing_interval/(24),tirz_av_15,'-.','Color',col_vec(4),'LineWidth',2);


ms=10;
lw=1.5;
plot([7 14],[tirz_av_5(1) tirz_av_5(8)],'o','Color',col_vec(2),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(2))
plot([7 14],[tirz_av_10(1) tirz_av_10(8)],'v','Color',col_vec(3),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(3))
plot([7 14],[tirz_av_15(1) tirz_av_15(8)],'diamond','Color',col_vec(4),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(4))


legend([p1,p2,p3],...
    "Tirzepatide 5 mg",...
    "Tirzepatide 10 mg",...
    "Tirzepatide 15 mg",...
    "Location","Northwest")

% legend([p1,p2,p3,p4,p5],...
%     "Tirzepatide 5 mg",...
%     "Tirzepatide 10 mg",...
%     "Tirzepatide 15 mg",...
%     "15\% weight loss",...
%     "20\% weight loss",...
%     "Location","Northwest")

xlabel("Dosing Interval (Days)")
ylabel("Percent Change in Body Weight")
grid on 
grid minor
% xlim([4,29])
xlim([6.5,28.5])
ylim([-25,0])
%ylabel("$\%\Delta$ in BW ($\%$)")
% title('Tirzepatide')
legend boxoff

xticks([7,10,14,17,21,24,28])
% xticklabels(['7','14','10','17','21','24','28'])

h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.XAxis.MinorTickValues = 7:1:28; % Minor ticks which don't line up with majors

%% panel B
subplot(2,2,2)
hold all

cost=dosing_interval/(24);
cost=cost(1)./cost;

plot(cost,tirz_av_5./tirz_av_5(1),'--','Color',col_vec(2),'LineWidth',2);
plot(cost,tirz_av_10./tirz_av_5(1),':','Color',col_vec(3),'LineWidth',2);
plot(cost,tirz_av_15./tirz_av_5(1),'-.','Color',col_vec(4),'LineWidth',2);
plot(cost,cost,'k','LineWidth',1)

ms=10;
lw=1.5;
plot([.5],[tirz_av_5(8)]./tirz_av_5(1),'o','Color',col_vec(2),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(2))
plot([.5],[tirz_av_10(8)]./tirz_av_5(1),'v','Color',col_vec(3),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(3))
plot([.5],[tirz_av_15(8)]./tirz_av_5(1),'diamond','Color',col_vec(4),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(4))
% plot([.5],[.5],'*k','MarkerSize',ms,'LineWidth',lw)

text(.41, tirz_av_5(8)./tirz_av_5(1), sprintf('%.0f\\%%', 100*tirz_av_5(8)./tirz_av_5(1)), 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');
text(.41, tirz_av_10(8)./tirz_av_5(1), sprintf('%.0f\\%%', 100*tirz_av_10(8)./tirz_av_5(1)), 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');
text(.39, tirz_av_15(8)./tirz_av_5(1), sprintf('%.0f\\%%', 100*tirz_av_15(8)./tirz_av_5(1)), 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');

legend(...
    "5 mg $\to$ 5 mg",...
    "5 mg $\to$ 10 mg",...
    "5 mg $\to$ 15 mg",...
    "Linear",...
    "Location","Southeast")

xlabel("Cost (rel.\ to 5 mg q1wk)")
ylabel("Efficacy (rel.\ to 5 mg q1wk)")
grid on 
grid minor
xlim([.25,1])
ylim([.25,1.5])
%ylabel("$\%\Delta$ in BW ($\%$)")
legend boxoff

xticks([.25,.5,.7,1])
yticks([.25,.5,.7,1,1.5])

h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
temp=flip(7./[7:1:28]);
h.XAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors
h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.YAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors

%% panel C
subplot(2,2,3)
hold all

cost=dosing_interval/(24);
cost=cost(1)./cost;

plot(cost,tirz_av_10./tirz_av_10(1),':','Color',col_vec(3),'LineWidth',2);
plot(cost,tirz_av_15./tirz_av_10(1),'-.','Color',col_vec(4),'LineWidth',2);
plot(cost,cost,'k','LineWidth',1)

ms=10;
lw=1.5;
plot([.5],[tirz_av_10(8)]./tirz_av_10(1),'v','Color',col_vec(3),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(3))
plot([.5],[tirz_av_15(8)]./tirz_av_10(1),'diamond','Color',col_vec(4),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(4))
% plot([.5],[.5],'*k','MarkerSize',ms,'LineWidth',lw)

text(.41, tirz_av_10(8)./tirz_av_10(1), sprintf('%.0f\\%%', 100*tirz_av_10(8)./tirz_av_10(1)), 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');
text(.41, tirz_av_15(8)./tirz_av_10(1), sprintf('%.0f\\%%', 100*tirz_av_15(8)./tirz_av_10(1)), 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');

legend(...
    "10 mg $\to$ 10 mg",...
    "10 mg $\to$ 15 mg",...
    "Linear",...
    "Location","Southeast")

xlabel("Cost (rel.\ to 10 mg q1wk)")
ylabel("Efficacy (rel.\ to 10 mg q1wk)")
grid on 
grid minor
xlim([.25,1])
ylim([.25,1.15])
%ylabel("$\%\Delta$ in BW ($\%$)")
legend boxoff

xticks([.25,.5,.7,1])
yticks([.25,.5,.7,1,1.15])

h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
temp=flip(7./[7:1:28]);
h.XAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors
h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.YAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors

%% panel D
subplot(2,2,4)
hold all

cost=dosing_interval/(24);
cost=cost(1)./cost;

plot(cost,tirz_av_15./tirz_av_15(1),'-.','Color',col_vec(4),'LineWidth',2);
plot(cost,cost,'k','LineWidth',1)

ms=10;
lw=1.5;
plot([.5],[tirz_av_15(8)]./tirz_av_15(1),'diamond','Color',col_vec(4),'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',col_vec(4))
% plot([.5],[.5],'*k','MarkerSize',ms,'LineWidth',lw)

text(.41, tirz_av_15(8)./tirz_av_15(1), sprintf('%.0f\\%%', 100*tirz_av_15(8)./tirz_av_15(1)), 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');

legend(...
    "15 mg $\to$ 15 mg",...
    "Linear",...
    "Location","Southeast")

xlabel("Cost (rel.\ to 15 mg q1wk)")
ylabel("Efficacy (rel.\ to 15 mg q1wk)")
grid on 
grid minor
xlim([.25,1])
ylim([.25,1])
%ylabel("$\%\Delta$ in BW ($\%$)")
legend boxoff

xticks([.25,.5,.7,1])
yticks([.25,.5,.7,1])

h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
temp=flip(7./[7:1:28]);
h.XAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors
h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.YAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors

%% format and save figure

fname="figTirz"

picturewidth = 40; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

% set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')

text(-.9, 2.05, 'A', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(0.12, 2.05, 'B', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(-.9, 1, 'C', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(0.12, 1, 'D', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');

legend boxoff
%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(hfig,strcat(fname,'.fig'))