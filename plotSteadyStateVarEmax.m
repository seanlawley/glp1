%Anil Cengiz - Nov 25,2024

%This scripts computes relative efficacies by taking weekly dosing as
%reference and using the data in semag_ss.mat

clear all
clc
close all

load("semag_ss.mat")

col_vec=["#e41a1c","#377eb8","#4daf4a","#984ea3"];

ref7days=holder(:,:,3);
holder_adj=holder./ref7days*100;

linwid=1.5;
hfig=figure()
hold on;

dosing_interval=dosing_interval(3:end);
holder_adj=holder_adj(:,:,3:end);


X=dosing_interval/(24);
X=X(1)./X;

for i=1:length(emaxi_vec)
    for  k=1:length(emaxs_vec)


        % scatter(dosing_interval/(24),squeeze(holder_adj(i,k,:)),80,'MarkerEdgeColor',col_vec(2),'MarkerFaceColor',col_vec(2),'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);
        p1=plot(X,squeeze(holder_adj(i,k,:))./100,'LineWidth',linwid);
        % alpha(.01)
    end
end

p2=plot(X,1.*X,'k--','LineWidth',linwid);
% legend([p1,p2],"Semaglutide 2.4 mg","Relative Frequency", "Location","Southeast")

legend([p1,p2],"Semaglutide 2.4 mg",...
    "Linear",...
    "Location","Southeast")

% text(.5, .75+.05,  '75\%', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');
% text(.5, .71-.01,  '71\%', 'FontSize', 35, 'Color', 'k', 'HorizontalAlignment', 'left');

xlabel("Cost (relative to once-weekly)")
ylabel("Efficacy (relative to once-weekly)")
grid on 
grid minor
xlim([.25,1])
ylim([.25,1])
xticks([.25,.5,.7,1])
yticks([.25,.5,.7,1])
h = gca; % Get axis to modify
h.XAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
temp=flip(7./[7:1:28]);
h.XAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors
h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.YAxis.MinorTickValues = temp; % Minor ticks which don't line up with majors



fname="figSteadyStateVarEmax"

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
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


