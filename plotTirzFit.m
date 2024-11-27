%% load data
close all
clear all
load('data_TirzFit.mat')

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

legend("Placebo","Tirzepatide 5 mg","Tirzepatide 10 mg","Tirzepatide 15 mg",'Location','southwest')
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