%% load data
close all
clear all

%% parameters
ps=0:.01:.5; % percentage of obese US adults

wl1=.17; % proportion body weight loss for q1wk
wl2=.12; % proportion body weight loss for q2wk
ys1=zeros(4,length(ps)); % hold calculations
ys2=ys1;
for i=1:length(ps)
    prob=ps(i);
    ys1(:,i)=calcBMI(prob,wl1);
    ys2(:,i)=calcBMI(2*prob,wl2);
end

S1_30=1e2.*ys1(1,:);
S1_35=1e2.*ys1(2,:);
S1_40=1e2.*ys1(3,:);
S2_30=1e2.*ys2(1,:);
S2_35=1e2.*ys2(2,:);
S2_40=1e2.*ys2(3,:);
livesSaved1=ys1(4,:);
livesSaved2=ys2(4,:);

dataTable=[...
    S1_30(26) S1_35(26) S1_40(26) livesSaved1(26);...
    S2_30(26) S2_35(26) S2_40(26) livesSaved2(26);...
    S1_30(51) S1_35(51) S1_40(51) livesSaved1(51);...
    S2_30(51) S2_35(51) S2_40(51) livesSaved2(51)];

round(dataTable(:,1:end-1))
dataTable(:,end)

%% convert ps to cost

numObeseAdults=0.41*260004955; % number of obese adults in US
mg=ps.*numObeseAdults*52*2.4./1e9;

%% Graphing

linwid=3;
ms=10;
hfig=figure();
hold on;

col_vec=["#1b9e77","#d95f02","#7570b3"];

subplot(2,2,1)
hold all
plot(mg,S1_30,'-','Color',col_vec(1),'LineWidth',linwid)
plot(mg,S2_30,':','Color',col_vec(2),'LineWidth',linwid)
plot([mg(26) mg(51)],[S1_30(26) S1_30(51)],'o','MarkerFaceColor',col_vec(1),'MarkerSize',ms)
plot([mg(26) mg(51)],[S2_30(26) S2_30(51)],'s','MarkerFaceColor',col_vec(2),'MarkerSize',ms)
xlabel("mg of semaglutide (billions/year)")
ylabel("\% with BMI$\ge30$")
legend("q1wk", "q2wk", "Location","Northeast")
legend boxoff
xlim([0,mg(end)])
grid on 
grid minor

subplot(2,2,2)
hold all
plot(mg,S1_35,'-','Color',col_vec(1),'LineWidth',linwid)
plot(mg,S2_35,':','Color',col_vec(2),'LineWidth',linwid)
plot([mg(26) mg(51)],[S1_35(26) S1_35(51)],'o','MarkerFaceColor',col_vec(1),'MarkerSize',ms)
plot([mg(26) mg(51)],[S2_35(26) S2_35(51)],'s','MarkerFaceColor',col_vec(2),'MarkerSize',ms)
xlabel("mg of semaglutide (billions/year)")
ylabel("\% with BMI$\ge35$")
legend("q1wk", "q2wk", "Location","Northeast")
legend boxoff
xlim([0,mg(end)])
grid on 
grid minor

subplot(2,2,3)
hold all
plot(mg,S1_40,'-','Color',col_vec(1),'LineWidth',linwid)
plot(mg,S2_40,':','Color',col_vec(2),'LineWidth',linwid)
plot([mg(26) mg(51)],[S1_40(26) S1_40(51)],'o','MarkerFaceColor',col_vec(1),'MarkerSize',ms)
plot([mg(26) mg(51)],[S2_40(26) S2_40(51)],'s','MarkerFaceColor',col_vec(2),'MarkerSize',ms)
xlabel("mg of semaglutide (billions/year)")
ylabel("\% with BMI$\ge40$")
legend("q1wk", "q2wk", "Location","Northeast")
legend boxoff
xlim([0,mg(end)])
grid on 
grid minor

subplot(2,2,4)
hold all
plot(mg,livesSaved1./1e3,'-','Color',col_vec(1),'LineWidth',linwid)
plot(mg,livesSaved2./1e3,':','Color',col_vec(2),'LineWidth',linwid)
plot([mg(26) mg(51)],[livesSaved1(26) livesSaved1(51)]./1e3,'o','MarkerFaceColor',col_vec(1),'MarkerSize',ms)
plot([mg(26) mg(51)],[livesSaved2(26) livesSaved2(51)]./1e3,'s','MarkerFaceColor',col_vec(2),'MarkerSize',ms)
xlabel("mg of semaglutide (billions/year)")
ylabel("Lives Saved ($10^3$/year)")
legend("q1wk", "q2wk", "Location","Southeast")
legend boxoff
xlim([0,mg(end)])
grid on 
grid minor

fname="figBMI"

picturewidth = 40; % set this parameter and keep it forever
hw_ratio = .8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

% set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[30 30 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')


text(-3.8*2.4*52/48, 600, 'A', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(-.5*2.4*52/48, 600, 'B', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(-3.8*2.4*52/48, 260, 'C', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
text(-.5*2.4*52/48, 260, 'D', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontName', 'Arial');

%print(hfig,fname,'-dpng','-painters')
%set(hfig, 'Position', get(0, 'Screensize'));
exportgraphics(hfig, strcat(fname,'.pdf'), 'ContentType', 'vector');
saveas(hfig,strcat(fname,'.fig'))