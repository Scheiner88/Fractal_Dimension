clear; clc; close all;
% -------------------------------------------------------------------------
load("InputData/Chembo_OutputData_beta_15-0_Fractal_2022_12_25-19_35_06.mat");

set_smdt = 1;

log_C_min = -1;
log_C_max = -0.5;
log_l_min = -10;
% -------------------------------------------------------------------------
log_C_comp = smoothdata(log_C,1,'gaussian',set_smdt);

log_C_comp(log_C < log_C_min) = NaN;
log_C_comp(log_C > log_C_max) = NaN;
log_C_comp(log_l < log_l_min,:) = NaN;

dlogC = diff(log_C_comp);
dlogl = diff(log_l);

difflogC = dlogC./dlogl;
difflogC(difflogC == inf) = NaN;

mean_difflogC = max(difflogC);

F = figure;
F.WindowState = 'maximized';
tiledlayout(2,2);

nexttile
plot(log_l,log_C,'LineWidth',2);
xlabel('log({\itr})'); ylabel('log({\itC})');
graph_setup(14);

nexttile
plot(log_l,log_C_comp,'LineWidth',2);
xlabel('log({\itr})'); ylabel('log({\itC})');
graph_setup(14);

nexttile
plot(m,mean_difflogC,'-','Marker','.','MarkerSize',15);
hold on;
plot(0:50,0:50,'--');
xlabel('{\itm}'); ylabel('{\itD}_2');
graph_setup(14);

nexttile
plot(log_l(1:end-1),difflogC(:,1:end),'-','LineWidth',2);
xlabel('log({\itr})'); ylabel('dlog({\itC})/dlog({\itr})');
graph_setup(14);

% figure
% plot(log_l,log_C,'LineWidth',2);
% xlabel('lg({\itr})'); ylabel('lg({\itC})');
% legend('{\itm} = ' + string(m),Location="southeast");
% graph_setup(14);
% 
% figure
% plot(log_l(1:end-1),difflogC(:,1:end),'-','LineWidth',2);
% xlabel('lg({\itr})'); ylabel('D^{1}[lg({\itC})]');
% legend('{\itm} = ' + string(m),Location="northeast");
% graph_setup(14);
% 
% D_2 = [0
% 4.644908953
% 9.226614049
% 12.51104539
% 15.0025795
% 16.48210802
% 18.16981165
% 19.85898532
% 21.42626765
% 22.99217387
% 24.3900713
% ];
% 
% figure
% plot([0 m],D_2,'-','Marker','.','MarkerSize',20,'LineWidth',2);
% y_lim_pl = ylim;
% hold on;
% plot(0:50,0:50,'--');
% ylim(y_lim_pl);
% xlabel('{\itm}'); ylabel('{\itD}_2');
% graph_setup(14);

D_2 = mean_difflogC';

function graph_setup(FontSize)
ab=findobj(gcf);
alltext=findall(ab,'Type','text');
allaxes=findall(ab,'Type','axes');
set(alltext,'FontName','Times New Roman', ...
    'FontWeight','Norm', ...
    'FontSize',FontSize);
set(allaxes,'FontName','Times New Roman', ...
    'FontWeight','Norm', ...
    'LineWidth',1,'FontSize',FontSize);
grid on; grid minor; box on;
end