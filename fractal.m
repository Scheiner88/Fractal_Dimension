clear; clc;
%-----------------------------НАСТРОЙКА------------------------------------
filename = 'InputData\Chembo_OutputData_beta_15-0';
save_output = 0; progres_bar = 1;

log_l = linspace(0.5,-1,500)'; m = 5:5:50;
max_length_S_fract = 1000; razr = 1;
work_tau_len = 5e6;
%---------------------------КОНЕЦ НАСТРОЙКИ--------------------------------
fwb = waitbar(0,'Loading your data');

load([filename,'.mat']);
date_str = char(datetime('now'),'yyyy_MM_dd-HH_mm_SS');

waitbar(0,fwb,'Processing your data');

S_lgc = false(max_length_S_fract,max(m),length(log_l));
S_lgc_all = zeros(1,length(m),length(log_l));
l_sq_3d = zeros(1,1,length(log_l)); l_sq_3d(1,1,:) = (10.^log_l).^2;
i_doub = round(max_length_S_fract/100); i_doub_proc = i_doub;

tau_idx = tau_mean_power(V,work_tau_len);

S = make_S(V,tau_idx,max(m),max_length_S_fract,razr);

tic

for i = 1:max_length_S_fract - 1
    S_cumsum = cumsum((S(i + 1:end,:) - S(i,:)).^2,2);
    S_lgc = lt(S_cumsum(:,m),l_sq_3d);
    S_lgc_all = S_lgc_all + sum(S_lgc,1);
    if i == i_doub_proc && progres_bar == 1
        i_fwb = i/max_length_S_fract;
        waitbar(i_fwb,fwb, ...
            ['Processing your data (',num2str(i_fwb*100),'%)']);
        i_doub_proc = i_doub_proc + i_doub;
    end
end
S_lgc_sum(:,:) = S_lgc_all(1,:,:); S_lgc_sum = (S_lgc_sum)';
log_C = log10(2*S_lgc_sum./(max_length_S_fract^2 - max_length_S_fract));

toc

close(fwb);

if save_output == 1
    save([filename,'_Fractal_',date_str],'m','max_length_S_fract', ...
        'tau_idx','razr','log_C','log_l');
end

figure
plot(log_l,log_C,'-','LineWidth',2);
xlabel('log({\itr})'); ylabel('log({\itC})');
legend('{\itm} = ' + string(m),'Location','southeastoutside');
graph_setup(14);

function tau_idx = tau_mean_power(x,work_tau_len)
x = x(end - work_tau_len:end);
n = length(x);
x_mean = mean(x(1:end - 1));
M = n - 1;
B_check = (1/M)*(sum((x(1:end - 1) - x_mean).*(x(2:end) - x_mean)));
for i = 2:n
    x_mean = mean(x(1:end - i));
    M = n - i;
    B = (1/M)*(sum((x(1:end - i) - x_mean).*(x(i + 1:end) - x_mean)));
    if B <= B_check/exp(1)
        tau_idx = i;
        break;
    end
end
end

function S = make_S(V,tau_idx,max_m,max_length_S_fract,razr)
x = fliplr(V(1:2*(max_length_S_fract*razr + (max_m - 1)*tau_idx)));
x = x - (1/2)*(max(x) + min(x)); x = (x/max(x) + 1)/2;
S(:,1) = x(1:end - (max_m - 1)*tau_idx);
S = [S(:,1),zeros(length(S(:,1)),max_m - 1)];
for i = 1:max_m - 1
    S(:,max_m - i + 1) = x(1 + (max_m - i)*tau_idx:end - (i - 1)*tau_idx);
end
S = S(1:razr:end,:); S = S(1:max_length_S_fract,:);
end

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