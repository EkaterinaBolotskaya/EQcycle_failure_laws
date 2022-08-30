%% Process and plot data for several failure laws

% 03/23/2022
% Ekaterina Bolotskaya

% This script plots earthquake cycle data for several failure laws. Data
% should be saved in Data_FL directory in .mat format. Run individual
% failure law scripts and save the data before running this one.

clearvars;
close all;

load ('Data_FL\DSWIS_sg.mat');
load ('Data_FL\DSWIS_gs.mat');
load ('Data_FL\SW.mat');
load ('Data_FL\RS.mat');
load ('Data_FL\EXP.mat');
load ('Data_FL\PAR.mat');

RS_full = RS;
id_s    = find(RS(:,2) == min(RS(:,2)),1);
RS      = RS(id_s:end,:);
RS(:,1) = RS(:,1)-RS(1,1);

%% Model parameters
Sn      = 20e6;                             % normal stress at depth of interest
K_s     = 1e6;                              % spring stiffness - factor < 1 is unstable
M       = 60e6;                             % mass

L       = 0.09;
a       = 0.0144;
b       = 0.02312;
A       = a*Sn;
B       = b*Sn;

V_0     = 1e-3;
V_pc    = 5*V_0;
v_st    = min(RS(:,2));
v_stc   = 400*v_st;

mu_s    = 0.7;
mu_d    = 0.6; 

tau_d   = Sn*mu_d;
d_tau   = Sn*(mu_s-mu_d);

mu_i    = 0.666519587834406;
tau_th  = mu_i*Sn;
mu_d_rs = 0.64395;

h       = -K_s/Sn+b/L;

% plotting parameters
lw      = 0.75;
fs      = 11;

V_int   = V_pc;                       % cross velocity plots at this value
D_int   = 1.8e-1;                     % cross slip plots at this value

% Colors
b_col   = [0.231372549019608 0.298039215686275 0.752941176470588];
r_col   = [0.705882352941177 0.015686274509804 0.149019607843137];
g_col   = [0 .7 0];

%% Allocate
ax100   = gobjects(2, 3);
ax101   = gobjects(2, 3);
ax102   = gobjects(2, 3);

%% Nondimentional variables for plotting
nd_t    = sqrt(M/K_s);
nd_u    = Sn*(mu_s-mu_d)/K_s;
nd_v    = nd_u/nd_t;
nd_a    = nd_v/nd_t;

%% Reference solution for static -> dynamic
% Full
N       = 5e3;
E_c     = Sn*(mu_s-mu_d)/M;
B_c     = K_s*V_0/M;
C_c     = K_s/M;
t_sd    = linspace(0, pi/sqrt(C_c), N);
u_sd    = E_c/C_c-E_c*cos(sqrt(C_c)*t_sd)/C_c-B_c*sin(sqrt(C_c)*t_sd)/(C_c^(3/2))+B_c/C_c*t_sd;
v_sd    = E_c*sin(sqrt(C_c)*t_sd)/sqrt(C_c)-B_c*cos(sqrt(C_c)*t_sd)/C_c+B_c/C_c;

%% Analytical solutions for RS 
% Shear rate = 0
t_EQ     = a/h*(1/v_stc);
time_a   = fliplr(t_EQ-(logspace(-10,log10(t_EQ),N)));
Sl       = -a/h*log(1-v_stc*h*time_a/a);           % slip
Sl_r     = (1/v_stc-h*time_a/a).^(-1);             % slip rate

% Shear rate = const ~= 0
tau_r    = K_s*V_0;                                % shear rate
t_EQ1    = a*Sn/tau_r*log(tau_r/h/Sn/v_st+1);      % time to instability  
time_a1  = fliplr(t_EQ1-(logspace(-10,log10(t_EQ1),N)));
Sl1      = -a/h*log(v_st*h*Sn/tau_r*(1-exp(tau_r*time_a1/a/Sn))+1);              % slip
Sl_r1    = ((1/v_st+h*Sn/tau_r)*exp(-tau_r*time_a1/A)-h*Sn/tau_r).^(-1);         % slip rate

%% find max V ind
im_RS          = find(RS(:,2)       == max(RS(:,2)),1);
im_RSa         = find(Sl_r1         == max(Sl_r1),1);
im_RSac        = find(Sl_r          == max(Sl_r),1);
im_EXP         = find(EXP(:,2)      == max(EXP(:,2)),1);
im_SW          = find(SW(:,2)       == max(SW(:,2)),1);
im_DSWIS_sg    = find(DSWIS_sg(:,2) == max(DSWIS_sg(:,2)),1);
im_C           = find(v_sd          == max(v_sd),1);
im_DSWIS_gs    = find(DSWIS_gs(:,2) == max(DSWIS_gs(:,2)),1);
im_PAR         = find(PAR(:,2)      == max(PAR(:,2)),1);

V_ei           = max([min(EXP(1:im_EXP,2)) min(PAR(1:im_PAR,2)) min(SW(1:im_SW,2)) min(DSWIS_sg(1:im_DSWIS_sg,2)) min(DSWIS_gs(1:im_DSWIS_gs,2))]);

%% find V crossing ind
[~,ic_RSa]     = min(abs(Sl_r1-V_int));
[~,ic_RSac]    = min(abs(Sl_r-V_int));
[~,ic_RS]      = min(abs(RS(1:im_RS,2)-V_int));
[~,ic_EXP]     = min(abs(EXP(1:im_EXP,2)-V_int));
[~,ic_SW]      = min(abs(SW(1:im_SW,2)-V_int));
ic_DSWIS_sg_ar = find(islocalmin(abs(DSWIS_sg(1:im_DSWIS_sg,2)-V_int)));
ic_DSWIS_sg    = ic_DSWIS_sg_ar(end);

ic_DSWIS_gs_ar = find(islocalmin(abs(DSWIS_gs(1:im_DSWIS_gs,2)-V_int)));
ic_DSWIS_gs    = ic_DSWIS_gs_ar(end);

ic_PAR_ar      = find(islocalmin(abs(PAR(1:im_PAR,2)-V_int)));
ic_PAR         = ic_PAR_ar(end);

[~,ic_C]       = min(abs(v_sd(1:im_C)-V_int));

%% Velocity vs. time reverse crossing at V_th
t_c            = RS(ic_RS,1);
t_m            = RS(im_RS,1);

t_RSa          = time_a1-time_a1(ic_RSa)+t_c-t_m;
t_RSac         = time_a-time_a(ic_RSac)+t_c-t_m;
t_RS           = RS(1:im_RS,1)-RS(ic_RS,1)+t_c-t_m;
t_EXP          = EXP(1:im_EXP,1)-EXP(ic_EXP,1)+t_c-t_m;
t_SW           = SW(1:im_SW,1)-SW(ic_SW,1)+t_c-t_m;
t_DSWIS_sg     = DSWIS_sg(1:im_DSWIS_sg,1)-DSWIS_sg(ic_DSWIS_sg,1)+t_c-t_m;
t_DSWIS_gs     = DSWIS_gs(1:im_DSWIS_gs,1)-DSWIS_gs(ic_DSWIS_gs,1)+t_c-t_m;
t_PAR          = PAR(1:im_PAR,1)-PAR(ic_PAR,1)+t_c-t_m;
t_C            = t_sd(1:im_C)-t_sd(ic_C)+t_c-t_m;

v_RSa          = Sl_r1;
v_RSac         = Sl_r;
v_RS           = RS(1:im_RS,2);
v_EXP          = EXP(1:im_EXP,2);
v_SW           = SW(1:im_SW,2);
v_DSWIS_sg     = DSWIS_sg(1:im_DSWIS_sg,2);
v_DSWIS_gs     = DSWIS_gs(1:im_DSWIS_gs,2);
v_PAR          = PAR(1:im_PAR,2);
v_C            = v_sd(1:im_C);

d_RSa          = Sl1;
d_RSac         = Sl;
d_RS           = RS(1:im_RS,3);
d_EXP          = EXP(1:im_EXP,3);
d_SW           = SW(1:im_SW,3);
d_DSWIS_sg     = DSWIS_sg(1:im_DSWIS_sg,3);
d_DSWIS_gs     = DSWIS_gs(1:im_DSWIS_gs,3);
d_PAR          = PAR(1:im_PAR,3);
d_C            = u_sd(1:im_C);

%% find D crossing ind
[~,icD_RSa]      = min(abs(Sl1-D_int));
[~,icD_RSac]     = min(abs(Sl-D_int));
[~,icD_RS]       = min(abs(RS(1:im_RS,3)-D_int));
[~,icD_EXP]      = min(abs(EXP(1:im_EXP,3)-D_int));
[~,icD_SW]       = min(abs(SW(1:im_SW,3)-D_int));
[~,icD_DSWIS_sg] = min(abs(DSWIS_sg(1:im_DSWIS_sg,3)-D_int));
[~,icD_DSWIS_gs] = min(abs(DSWIS_gs(1:im_DSWIS_gs,3)-D_int));
[~,icD_PAR]      = min(abs(PAR(1:im_PAR,3)-D_int));
[~,icD_C]        = min(abs(u_sd(1:im_C)-D_int));

%% Slip vs. time reverse crossing at D_int
t_cD             = DSWIS_gs(icD_DSWIS_gs,1);
t_mD             = DSWIS_gs(im_DSWIS_gs,1);

tD_RSa           = time_a1-time_a1(icD_RSa)+t_cD-t_mD;
tD_RSac          = time_a-time_a(im_RSac);
tD_RS            = RS(1:im_RS,1)-RS(icD_RS,1)+t_cD-t_mD;
tD_EXP           = EXP(1:im_EXP,1)-EXP(icD_EXP,1)+t_cD-t_mD;
tD_SW            = SW(1:im_SW,1)-SW(icD_SW,1)+t_cD-t_mD;
tD_DSWIS_sg      = DSWIS_sg(1:im_DSWIS_sg,1)-DSWIS_sg(icD_DSWIS_sg,1)+t_cD-t_mD;
tD_DSWIS_gs      = DSWIS_gs(1:im_DSWIS_gs,1)-DSWIS_gs(icD_DSWIS_gs,1)+t_cD-t_mD;
tD_PAR           = PAR(1:im_PAR,1)-PAR(icD_PAR,1)+t_cD-t_mD;
tD_C             = t_sd(1:im_C)-t_sd(icD_C)+t_cD-t_mD;

%% Phase diagrams
figure(100);
set(gca,'FontSize',fs-2);
for i = 1:6
    if i == 1
        ax100(i) = subplot(2,3,i);
        hold on;
        set(gca, 'XScale', 'log');
        box on;
        title('RS', 'Interpreter', 'latex', 'FontSize', fs+1);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = RS_full;
        
        id_ei = find(var(2:end,4)>tau_th, 1, 'first');
        id_ip = find(var(2:end,2)>V_0, 1, 'first');
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
        id_ei2 = find(var(:,2)>V_ei, 1, 'last');
    elseif i == 2 
        ax100(i) = subplot(2,3,i);
        hold on;
        set(gca, 'XScale', 'log');
        box on;
        title('SW', 'Interpreter', 'latex', 'FontSize', fs+1);
        var = SW;
        im_var = im_SW;

        id_ei = find(var(2:end,2)>V_ei, 1, 'first');
        id_ip_ar = find(islocalmin(abs(var(1:im_var,2)-V_0)));
        id_ip    = id_ip_ar(end);
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
     elseif i == 3 
        ax100(i) = subplot(2,3,i);
        hold on;
        set(gca, 'XScale', 'log');
        box on;
        title('ECZ', 'Interpreter', 'latex', 'FontSize', fs+1);
        var = EXP;
        im_var = im_EXP;

        id_ei = find(var(2:end,2)>V_ei, 1, 'first');
        id_ip_ar = find(islocalmin(abs(var(1:im_var,2)-V_0)));
        id_ip    = id_ip_ar(end);
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
    elseif i == 4
        ax100(i) = subplot(2,3,i);
        hold on;
        set(gca, 'XScale', 'log');
        box on;
        title('DSWIS 1', 'Interpreter', 'latex', 'FontSize', fs+1);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = DSWIS_sg;
        im_var = im_DSWIS_sg;
        
        id_ei = find(var(2:end,2)>V_ei, 1, 'first');
        id_ip_ar = find(islocalmin(abs(var(1:im_var,2)-V_0)));
        id_ip    = id_ip_ar(end);
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
    elseif i == 5
        ax100(i) = subplot(2,3,i);
        hold on;
        set(gca, 'XScale', 'log');
        box on;
        title('DSWIS 2', 'Interpreter', 'latex', 'FontSize', fs+1);
        xlabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = DSWIS_gs;
        im_var = im_DSWIS_gs;
        
        id_ei = find(var(2:end,2)>V_ei, 1, 'first');
        id_ip_ar = find(islocalmin(abs(var(1:im_var,2)-V_0)));
        id_ip    = id_ip_ar(end);
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
    elseif i == 6
        ax100(i) = subplot(2,3,i);
        hold on;
        set(gca, 'XScale', 'log');
        box on;
        title('PCZ', 'Interpreter', 'latex', 'FontSize', fs+1);
        xlabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = PAR;
        im_var = im_PAR;
        
        id_ei = find(var(2:end,2)>V_ei, 1, 'first');
        id_ip_ar = find(islocalmin(abs(var(1:im_var,2)-V_0)));
        id_ip    = id_ip_ar(end);
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
    end

    semilogx((var(1:id_ei,2)/nd_v), (var(1:id_ei,4) - tau_d)/d_tau, 'k', 'LineWidth', lw);
    hold on;
    semilogx((var(id_ei:id_ip,2)/nd_v), (var(id_ei:id_ip,4) - tau_d)/d_tau, 'Color', r_col, 'LineWidth', lw);
    semilogx((var(id_ip:id_pc,2)/nd_v), (var(id_ip:id_pc,4) - tau_d)/d_tau, 'Color', g_col, 'LineWidth', lw);
    semilogx((var(id_pc:end,2)/nd_v), (var(id_pc:end,4) - tau_d)/d_tau, 'Color', b_col, 'LineWidth', lw);
    
    semilogx((var(1:id_ei,2)/nd_v), (var(1:id_ei,5) - tau_d)/d_tau, 'k--', 'LineWidth', lw);
    semilogx((var(id_ei:id_ip,2)/nd_v), (var(id_ei:id_ip,5) - tau_d)/d_tau, '--', 'Color', r_col, 'LineWidth', lw);
    semilogx((var(id_ip:id_pc,2)/nd_v), (var(id_ip:id_pc,5) - tau_d)/d_tau, '--', 'Color', g_col, 'LineWidth', lw);
    semilogx((var(id_pc:end,2)/nd_v), (var(id_pc:end,5) - tau_d)/d_tau, '--', 'Color', b_col, 'LineWidth', lw);
    
    semilogx((V_ei/nd_v)*ones(1,100), 1.1*(linspace(min((var(:,5) - tau_d)/d_tau), max((var(:,5) - tau_d)/d_tau), 100)), 'k-.', 'LineWidth', lw/2);
    semilogx((V_0/nd_v)*ones(1,100), 1.1*(linspace(min((var(:,5) - tau_d)/d_tau), max((var(:,5) - tau_d)/d_tau), 100)), 'k-.', 'LineWidth', lw/2);
    semilogx((V_pc/nd_v)*ones(1,100), 1.1*(linspace(min((var(:,5) - tau_d)/d_tau), max((var(:,5) - tau_d)/d_tau), 100)), 'k-.', 'LineWidth', lw/2);
    
    text(V_0*1.1/nd_v,0.07, '$V_0$', 'Interpreter', 'latex', 'FontSize', fs-3);
    text(V_pc*1.1/nd_v,0.07, '$V_{p-c}$', 'Interpreter', 'latex', 'FontSize', fs-3);
    text(V_ei*1.1/nd_v,0.07, '$V_{ei}$', 'Interpreter', 'latex', 'FontSize', fs-3);
    legend({'Early interseismic', 'Late interseismic', 'Preseismic', 'Coseismic'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
    grid on;
    xlim([0.85*min(var(:,2))/nd_v 1.15*max(var(:,2))/nd_v]);
    ylim([min((var(:,5) - tau_d)/d_tau) max((var(:,5) - tau_d)/d_tau)]*1.04);
    
end

linkaxes(ax100, 'xy');
hold off;

%% Energy curves
figure(101);
set(gca,'FontSize',fs-2);
for i = 1:6
    if i == 1
        ax101(i) = subplot(2,3,i);
        hold on;
        title('RS', 'Interpreter', 'latex', 'FontSize', fs+1);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = RS;
    
        idx_w = find(var(:,4)>tau_d, 1, 'first');
        u_int = interp1([var(idx_w-1,4) var(idx_w,4)], [var(idx_w-1,3) var(idx_w,3)], tau_d, 'linear');
        
        idm = find(var(:,4)==max(var(:,4)));
        idk = find(var(:,5)>=var(:,4), 1, 'last');
        idd = find(var(1:idk,4)>tau_d, 1, 'last');
    elseif i == 2 
        ax101(i) = subplot(2,3,i);
        hold on;
        title('SW', 'Interpreter', 'latex', 'FontSize', fs+1);
        var = SW;

        idx_w = find(var(:,4)>tau_d, 1, 'first');
        u_int = interp1([var(idx_w-1,4) var(idx_w,4)], [var(idx_w-1,3) var(idx_w,3)], tau_d, 'linear');
        
        idm = find(var(:,4)==max(var(:,4)));
        idk = find(var(:,5)>=tau_d, 1, 'last');
        idd = find(var(1:idk,4)>tau_d, 1, 'last');
     elseif i == 3 
        ax101(i) = subplot(2,3,i);
        hold on;
        title('ECZ', 'Interpreter', 'latex', 'FontSize', fs+1);
        var = EXP;

        idx_w = find(var(:,4)>tau_d, 1, 'first');
        u_int = interp1([var(idx_w-1,4) var(idx_w,4)], [var(idx_w-1,3) var(idx_w,3)], tau_d, 'linear');
        
        idm = find(var(:,4)==max(var(:,4)));
        idk = find(var(:,5)>=tau_d, 1, 'last');
        idd = find(var(1:idk,4)>tau_d, 1, 'last');
    elseif i == 4
        ax101(i) = subplot(2,3,i);
        hold on;
        title('DSWIS 1', 'Interpreter', 'latex', 'FontSize', fs+1);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = DSWIS_sg;

        idx_w = find(var(:,4)>tau_d, 1, 'first');
        u_int = interp1([var(idx_w-1,4) var(idx_w,4)], [var(idx_w-1,3) var(idx_w,3)], tau_d, 'linear');
        
        idm = find(var(:,4)==max(var(:,4)));
        idk = find(var(:,5)>=tau_d, 1, 'last');
        idd = find(var(1:idk,4)>tau_d, 1, 'last');
    elseif i == 5
        ax101(i) = subplot(2,3,i);
        hold on;
        title('DSWIS 2', 'Interpreter', 'latex', 'FontSize', fs+1);
        xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = DSWIS_gs;

        idx_w = find(var(:,4)>tau_d, 1, 'first');
        u_int = interp1([var(idx_w-1,4) var(idx_w,4)], [var(idx_w-1,3) var(idx_w,3)], tau_d, 'linear');
        
        idm = find(var(:,4)==max(var(:,4)));
        idk = find(var(:,5)>=tau_d, 1, 'last');
        idd = find(var(1:idk,4)>tau_d, 1, 'last');
    elseif i == 6
        ax101(i) = subplot(2,3,i);
        hold on;
        title('PCZ', 'Interpreter', 'latex', 'FontSize', fs+1);
        xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = PAR;

        idx_w = find(var(:,4)>tau_d, 1, 'first');
        u_int = interp1([var(idx_w-1,4) var(idx_w,4)], [var(idx_w-1,3) var(idx_w,3)], tau_d, 'linear');
        
        idm = find(var(:,4)==max(var(:,4)));
        idk = find(var(:,5)>=tau_d, 1, 'last');
        idd = find(var(1:idk,4)>tau_d, 1, 'last');
    end
       
    xr = [u_int; var(idx_w:idm,3)]/nd_u;
    y1r = ([mu_d*Sn; var(idx_w:idm,4)]-mu_d*Sn)/Sn/(mu_s-mu_d);
    y2r = zeros(size(y1r));
    fill([xr; flipud(xr)], [y1r; flipud(y2r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
    hold on;
    xb = var(idm:idd,3)/nd_u;
    y1b = (var(idm:idd,4)-mu_d*Sn)/Sn/(mu_s-mu_d);
    y2b = zeros(size(y1b));
    fill([xb; flipud(xb)], [y1b; flipud(y2b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
    xg         = var(idm:idk,3)/nd_u;
    y1g        = (var(idm:idk,5)-mu_d*Sn)/Sn/(mu_s-mu_d);
    y2g        = (var(idm:idk,4)-mu_d*Sn)/Sn/(mu_s-mu_d);
    p3 = fill([xg; flipud(xg)], [y2g; flipud(y1g)], 'k', 'EdgeColor', [.4 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.2);
       
    p1 = plot(var(:,3)/nd_u, (var(:,4)-tau_d)/d_tau, 'k', 'LineWidth', lw);
    p2 = plot(var(:,3)/nd_u, (var(:,5)-tau_d)/d_tau, 'k--', 'LineWidth', lw);
    legend([p1 p2], {'Friction', 'Loading'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northeast');
    
    xlim([0 2]);
    ylim([min((var(:,5) - tau_d)/d_tau) max((var(:,5) - tau_d)/d_tau)]*1.04);
    box on;
    grid on;

end

linkaxes(ax101, 'xy');
hold off;

%% Slip and slip rate single cycle
fig102 = figure(102);
set(gca,'FontSize',fs-2);
set(fig102, 'defaultAxesColorOrder',[r_col; b_col]);

for i = 1:6
    if i == 1
        ax102(i) = subplot(2,3,i);
        hold on;
        yyaxis right;
        set(gca, 'YScale', 'log');
        title('RS', 'Interpreter', 'latex', 'FontSize', fs+1);
        var = RS_full;
    elseif i == 2 
        ax102(i) = subplot(2,3,i);
        hold on;
        yyaxis right;
        set(gca, 'YScale', 'log');
        title('SW', 'Interpreter', 'latex', 'FontSize', fs+1);
        var = SW;
     elseif i == 3 
        ax102(i) = subplot(2,3,i);
        hold on;
        yyaxis right;
        set(gca, 'YScale', 'log');
        title('ECZ', 'Interpreter', 'latex', 'FontSize', fs+1);
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = EXP;
    elseif i == 4
        ax102(i) = subplot(2,3,i);
        hold on;
        yyaxis right;
        set(gca, 'YScale', 'log');
        title('DSWIS 1', 'Interpreter', 'latex', 'FontSize', fs+1);
        xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = DSWIS_sg;
    elseif i == 5
        ax102(i) = subplot(2,3,i);
        hold on;
        yyaxis right;
        set(gca, 'YScale', 'log');
        title('DSWIS 2', 'Interpreter', 'latex', 'FontSize', fs+1);
        xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = DSWIS_gs;
    elseif i == 6
        ax102(i) = subplot(2,3,i);
        hold on;
        yyaxis right;
        set(gca, 'YScale', 'log');
        title('PCZ', 'Interpreter', 'latex', 'FontSize', fs+1);
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        var = PAR;
    end
       
    yyaxis right;
    p(1) = semilogy(var(:,1)/nd_t, var(:,3)/nd_u, 'LineWidth', lw);
    hold on;
    p(2) = semilogy(var(:,1)/nd_t, V_0*var(:,1)/nd_u, 'LineWidth', lw);
    set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0]);
    hold off;
    
    yyaxis left;
    p(3) = semilogy(var(:,1)/nd_t, var(:,2)/nd_v, 'LineWidth', lw);
    if (i == 1) || (i == 4)
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    end
    hold on;
    semilogy(var(:,1)/nd_t, V_0/nd_v*ones(size(var(:,1))), 'k-.', 'LineWidth', lw/2);
    semilogy(var(:,1)/nd_t, V_pc/nd_v*ones(size(var(:,1))), 'k-.', 'LineWidth', lw/2);
    semilogy(var(:,1)/nd_t, V_ei/nd_v*ones(size(var(:,1))), 'k-.', 'LineWidth', lw/2);
    
    text(200,V_0*2/nd_v, '$V_0$', 'Interpreter', 'latex', 'FontSize', fs-2);
    text(200,V_pc*2/nd_v, '$V_{p-c}$', 'Interpreter', 'latex', 'FontSize', fs-2);
    text(200,V_ei*2/nd_v, '$V_{ei}$', 'Interpreter', 'latex', 'FontSize', fs-2);
    
    legend([p(1) p(2)],'Slider','Load point', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0]);

end

yyaxis left;
linkaxes(ax102, 'xy');
hold off;

%% Plots by phase
cyc_d   = zeros(1,7);
cum_u   = zeros(1,7);
tau_sta = zeros(1,7);
max_cv  = zeros(1,7);
ie_t    = zeros(1,7);
li_t    = zeros(1,7);
pr_t    = zeros(1,7);
co_t    = zeros(1,7);
ie_u    = zeros(1,7);
li_u    = zeros(1,7);
pr_u    = zeros(1,7);
co_u    = zeros(1,7);


for i = 1:7
    if i == 1 
        var = [t_sd; v_sd; u_sd]';
        col = [0 0 0];
        idpl = 7;

        id_ei = find(var(2:end,2)>1.0001*V_ei, 1, 'first');
        id_ip = find(var(2:end,2)>V_0, 1, 'first');
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
        id_fl = id_pc;

        wr_id = 7;
    elseif i == 2
        var = DSWIS_sg;
        col = [0.4660, 0.6740, 0.1880];
        idpl = 1;

        id_ei = find(var(2:end,2)>1.0001*V_ei, 1, 'first');
        id_ip = find(var(2:end,2)>V_0, 1, 'first');
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
        id_fl = id_pc+find(var(id_pc:end,4)>tau_d, 1, 'last');

        wr_id = 4;
    elseif i == 3
        var = DSWIS_gs;
        col = [0.3010, 0.7450, 0.9330];
        idpl = 2;

        id_ei = find(var(2:end,2)>1.0001*V_ei, 1, 'first');
        id_ip = find(var(2:end,2)>V_0, 1, 'first');
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
        id_fl = id_pc+find(var(id_pc:end,4)>tau_d, 1, 'last');

        wr_id = 5;
    elseif i == 4
        var = PAR;
        col = [0.6350, 0.0780, 0.1840];
        idpl = 3;

        id_ei = find(var(2:end,2)>1.0001*V_ei, 1, 'first');
        id_ip_ar = find(islocalmin(abs(var(1:im_var,2)-V_0)));
        id_ip    = id_ip_ar(end);
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
        id_fl = id_pc+find(var(id_pc:end,4)>tau_d, 1, 'last');

        wr_id = 6;
    elseif i == 5
        var = EXP;
        col = [0, 0.4470, 0.7410];
        idpl = 4;

        id_ei = find(var(2:end,2)>1.0001*V_ei, 1, 'first');
        id_ip = find(var(2:end,2)>V_0, 1, 'first');
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
        id_fl = id_pc+find(var(id_pc:end,4)>tau_d*1.001, 1, 'last');

        wr_id = 3;
    elseif i == 6
        var = RS_full;
        col = [0.8500, 0.3250, 0.0980];
        idpl = 5;

        id_ei = find(var(2:end,4)>tau_th, 1, 'first');
        id_ip = find(var(2:end,2)>V_0, 1, 'first');
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
        id_fl = length(var(:,2));

        wr_id = 1;
    elseif i == 7
        var = SW;
        col = [0.9290, 0.6940, 0.1250];
        idpl = 6;

        id_ei = find(var(2:end,2)>1.0001*V_ei, 1, 'first');
        id_ip = find(var(2:end,2)>V_0, 1, 'first');
        id_pc = find(var(:,2)>=V_pc, 1, 'first');
        id_fl = id_pc+find(var(id_pc:end,4)>tau_d, 1, 'last');

        wr_id = 2;
    end

    int_tip = interp1(var(id_ip-1:id_ip+1,2), var(id_ip-1:id_ip+1,1), V_0);
    if i == 4
        int_tip = var(id_ip,1);
        v_int   = var(id_ip,2);
    else
        v_int = V_0;
    end
    int_uip = interp1(var(id_ip-1:id_ip+1,2), var(id_ip-1:id_ip+1,3), V_0);
    
    int_tpc = interp1(var(id_pc-1:id_pc+1,2), var(id_pc-1:id_pc+1,1), V_pc);
    int_upc = interp1(var(id_pc-1:id_pc+1,2), var(id_pc-1:id_pc+1,3), V_pc);

    
    cum_u(wr_id)   = max(var(:,3));
    cyc_d(wr_id)   = max(var(:,1));
    if i~=1
        tau_sta(wr_id) = max(var(:,4))-min(var(:,4));
        
    else
        tau_sta(wr_id) = 2*Sn*(mu_s-mu_d);
        
    end
    max_cv(wr_id)  = max(var(:,2));
    ie_t(wr_id)    = var(id_ei,1)-var(1,1);
    li_t(wr_id)    = var(id_ip,1)-var(id_ei,1);
    pr_t(wr_id)    = var(id_pc,1)-var(id_ip,1);
    co_t(wr_id)    = var(end,1)-var(id_pc,1);
    ie_u(wr_id)    = var(id_ei,3)-var(1,3);
    li_u(wr_id)    = var(id_ip,3)-var(id_ei,3);
    pr_u(wr_id)    = var(id_pc,3)-var(id_ip,3);
    co_u(wr_id)    = var(end,3)-var(id_pc,3);
   
    % Late interseismic phase
    fig1111 = figure(1111);
    set(gca,'FontSize',fs-2);
    
    subplot(121)
    p1(idpl) = plot(([var(id_ei:id_ip,1); int_tip]-var(id_ei,1))/nd_t, [var(id_ei:id_ip,2); v_int]/nd_v, 'Color', col, 'LineWidth', lw);
    hold on;
    
    subplot(122)
    p1(idpl) = plot(([var(id_ei:id_ip,1); int_tip]-var(id_ei,1))/nd_t, ([var(id_ei:id_ip,3); int_uip]-var(id_ei,3))/nd_u, 'Color', col,'LineWidth', lw);
    hold on;

    % Preseismic phase
    fig1112 = figure(1112);
    set(gca,'FontSize',fs-2);
    subplot(121)
    p2(idpl) = plot(([int_tip; var(id_ip+1:id_pc-1,1); int_tpc]-int_tip)/nd_t, [v_int; var(id_ip+1:id_pc-1,2); V_pc]/nd_v, 'Color', col, 'LineWidth', lw);
    hold on;
    
    subplot(122)
    p2(idpl) = plot(([int_tip; var(id_ip+1:id_pc-1,1); int_tpc]-int_tip)/nd_t, ([int_uip; var(id_ip+1:id_pc-1,3); int_upc]-var(id_ip,3))/nd_u, 'Color', col, 'LineWidth', lw);
    hold on;
    
    % Coseismic phase
    fig1113 = figure(1113);
    set(gca,'FontSize',fs-2);
    subplot(121)
    p3(idpl) = plot((var(id_pc:id_fl,1)-var(id_pc,1))/nd_t, var(id_pc:id_fl,2)/nd_v, 'Color', col, 'LineWidth', 1.1);
    hold on;
    plot((var(id_fl:end,1)-var(id_pc,1))/nd_t, var(id_fl:end,2)/nd_v, '--','Color', col, 'LineWidth', lw);
    
    subplot(122)
    p3(idpl) = plot((var(id_pc:id_fl,1)-var(id_pc,1))/nd_t, (var(id_pc:id_fl,3)-var(id_pc,3))/nd_u, 'Color', col, 'LineWidth', 1.1);
    hold on;
    plot((var(id_fl:end,1)-var(id_pc,1))/nd_t, (var(id_fl:end,3)-var(id_pc,3))/nd_u, '--','Color', col, 'LineWidth', lw);
   

end

figure(1111)
subplot(121)
yline(V_0/nd_v, 'k-.', 'LineWidth', lw/2);
text(160,V_0*1.035/nd_v, '$V_0$', 'Interpreter', 'latex', 'FontSize', fs-2);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend([p1(5) p1(6) p1(4) p1(1) p1(2) p1(3) p1(7)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northeast');
grid on;
axis tight;

subplot(122)
ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
grid on;
axis tight;

sgtitle('Late interseimic phase', 'Interpreter', 'latex', 'FontSize', fs+1);

figure(1112)
subplot(121)
yline(V_pc/nd_v, 'k-.', 'LineWidth', lw/2);
text(5,V_pc*1.025/nd_v, '$V_{p-c}$', 'Interpreter', 'latex', 'FontSize', fs-2);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
grid on;
axis tight;
ylim([V_0/nd_v V_pc*1.05/nd_v]);

subplot(122)
ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend([p2(5) p2(6) p2(4) p2(1) p2(2) p2(3) p2(7)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northwest');
grid on;
axis tight;

sgtitle('Preseimic phase', 'Interpreter', 'latex', 'FontSize', fs+1);

figure(1113)
subplot(121)
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend([p3(5) p3(6) p3(4) p3(1) p3(2) p3(3) p3(7)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northwest');
grid on;
axis tight;

subplot(122)
ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
grid on;
axis tight;

sgtitle('Coseimic phase', 'Interpreter', 'latex', 'FontSize', fs+1);

%% Log slip rate vs. linear time - nondim
fig3 = figure(103);
set(gca,'FontSize',fs-2);
semilogy(abs(t_RSa(t_RSa<0))/nd_t, V_0*ones(size(abs(t_RSa(t_RSa<0))))/nd_v, 'k-.', 'LineWidth', lw/2);
hold on;
semilogy(abs(t_RSa(t_RSa<0))/nd_t, V_ei*ones(size(abs(t_RSa(t_RSa<0))))/nd_v, 'k-.', 'LineWidth', lw/2);
semilogy(abs(t_RSa(t_RSa<0))/nd_t, V_pc*ones(size(abs(t_RSa(t_RSa<0))))/nd_v, 'k-.', 'LineWidth', lw/2);
p(7) = semilogy(abs(t_C(t_C<0))/nd_t, v_C(t_C<0)/nd_v, 'k', 'LineWidth', lw);

p(1) = semilogy(abs(t_DSWIS_sg(t_DSWIS_sg<0))/nd_t, v_DSWIS_sg(t_DSWIS_sg<0)/nd_v, 'LineWidth', lw);
p(2) = semilogy(abs(t_DSWIS_gs(t_DSWIS_gs<0))/nd_t, v_DSWIS_gs(t_DSWIS_gs<0)/nd_v, 'LineWidth', lw);
p(3) = semilogy(abs(t_PAR(t_PAR<0))/nd_t, v_PAR(t_PAR<0)/nd_v, 'LineWidth', lw);
p(4) = semilogy(abs(t_EXP(t_EXP<0))/nd_t, v_EXP(t_EXP<0)/nd_v, 'LineWidth', lw);
p(5) = semilogy(abs(t_RS(t_RS<0))/nd_t, v_RS(t_RS<0)/nd_v, 'LineWidth', lw);
p(6) = semilogy(abs(t_SW(t_SW<0))/nd_t, v_SW(t_SW<0)/nd_v, 'LineWidth', lw);
text(140,V_0*1.5/nd_v, '$V_0$', 'Interpreter', 'latex', 'FontSize', fs-2);
text(140,V_pc*1.5/nd_v, '$V_{p-c}$', 'Interpreter', 'latex', 'FontSize', fs-2);
text(140,V_ei*1.5/nd_v, '$V_{ei}$', 'Interpreter', 'latex', 'FontSize', fs-2);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$ to instability', 'Interpreter', 'latex', 'FontSize', fs);
legend([p(5) p(6) p(4) p(1) p(2) p(3) p(7)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northeast');
grid on;
xlim([0.03/nd_t 0.55*RS(im_RS,1)/nd_t]);
ylim([0.99*v_st/nd_v 1.01*max(v_C)/nd_v]);

%% Log time vs. log vel - nondim
fig4 = figure(104);
set(gca,'FontSize',fs-2);
t_d = linspace(1e-9,10*max(abs(t_RSa(t_RSa<0)))/nd_t,1000);
loglog(t_d, V_0*ones(size(t_d))/nd_v, 'k-.', 'LineWidth', lw/2);
hold on;
loglog(t_d, V_ei*ones(size(t_d))/nd_v, 'k-.', 'LineWidth', lw/2);
loglog(t_d, V_pc*ones(size(t_d))/nd_v, 'k-.', 'LineWidth', lw/2);
p(7) = loglog(abs(t_C(t_C<0))/nd_t, v_C(t_C<0)/nd_v, 'k', 'LineWidth', lw);

p(1) = loglog(abs(t_DSWIS_sg(t_DSWIS_sg<0))/nd_t, v_DSWIS_sg(t_DSWIS_sg<0)/nd_v, 'LineWidth', lw);
p(2) = loglog(abs(t_DSWIS_gs(t_DSWIS_gs<0))/nd_t, v_DSWIS_gs(t_DSWIS_gs<0)/nd_v, 'LineWidth', lw);
p(3) = loglog(abs(t_PAR(t_PAR<0))/nd_t, v_PAR(t_PAR<0)/nd_v, 'LineWidth', lw);
p(4) = loglog(abs(t_EXP(t_EXP<0))/nd_t, v_EXP(t_EXP<0)/nd_v, 'LineWidth', lw);
p(5) = loglog(abs(t_RS(t_RS<0))/nd_t, v_RS(t_RS<0)/nd_v, 'LineWidth', lw);
p(6) = loglog(abs(t_SW(t_SW<0))/nd_t, v_SW(t_SW<0)/nd_v, 'LineWidth', lw);
text(2,V_0*1.5/nd_v, '$V_0$', 'Interpreter', 'latex', 'FontSize', fs-2);
text(2,V_pc*1.5/nd_v, '$V_{p-c}$', 'Interpreter', 'latex', 'FontSize', fs-2);
text(2,V_ei*1.5/nd_v, '$V_{ei}$', 'Interpreter', 'latex', 'FontSize', fs-2);

ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$ to instability', 'Interpreter', 'latex', 'FontSize', fs);
legend([p(5) p(6) p(4) p(1) p(2) p(3) p(7)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
grid on;
axis tight;
% xlim([2e-1 1.01*RS(im_RS,1)/nd_t]);
xlim([0.5e-1 1.01*RS(im_RS,1)/nd_t]);

%% Log time vs. log slip - nondim
fig5 = figure(105);
set(gca,'FontSize',fs-2);
t_d = linspace(1e-9,10*max(abs(t_RSa(t_RSa<0)))/nd_t,1000);
loglog([0 0], [0 0], 'k--', 'LineWidth', lw/2);
hold on;
plot([0 0], [0 0], 'k--', 'LineWidth', lw/2);
plot([0 0], [0 0], 'k--', 'LineWidth', lw/2);
p(7) = loglog(abs(tD_C(t_C<0))/nd_t, d_C(t_C<0)/nd_u, 'k', 'LineWidth', lw);

p(1) = loglog(abs(tD_DSWIS_sg(t_DSWIS_sg<0))/nd_t, d_DSWIS_sg(t_DSWIS_sg<0)/nd_u, 'LineWidth', lw);
p(2) = loglog(abs(tD_DSWIS_gs(t_DSWIS_gs<0))/nd_t, d_DSWIS_gs(t_DSWIS_gs<0)/nd_u, 'LineWidth', lw);
p(3) = loglog(abs(tD_PAR(t_PAR<0))/nd_t, d_PAR(t_PAR<0)/nd_u, 'LineWidth', lw);
p(4) = loglog(abs(tD_EXP(t_EXP<0))/nd_t, d_EXP(t_EXP<0)/nd_u, 'LineWidth', lw);
p(5) = loglog(abs(tD_RS(t_RS<0))/nd_t, d_RS(t_RS<0)/nd_u, 'LineWidth', lw);
p(6) = loglog(abs(tD_SW(t_SW<0))/nd_t, d_SW(t_SW<0)/nd_u, 'LineWidth', lw);

ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$ to instability', 'Interpreter', 'latex', 'FontSize', fs);
legend([p(5) p(6) p(4) p(1) p(2) p(3) p(7)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
grid on;
axis tight;
xlim([1 1.01*RS(im_RS,1)/nd_t]);
ylim([1e-5 1.3]);

%% Time vs. log slip - nondim
fig6 = figure(106);
set(gca,'FontSize',fs-2);
semilogy([0 0], [0 0], 'k--', 'LineWidth', lw/2);
hold on;
plot([0 0], [0 0], 'k--', 'LineWidth', lw/2);
plot([0 0], [0 0], 'k--', 'LineWidth', lw/2);
p(7) = semilogy(abs(tD_C(t_C<0))/nd_t, d_C(t_C<0)/nd_u, 'k', 'LineWidth', lw);

p(1) = semilogy(abs(tD_DSWIS_sg(t_DSWIS_sg<0))/nd_t, d_DSWIS_sg(t_DSWIS_sg<0)/nd_u, 'LineWidth', lw);
p(2) = semilogy(abs(tD_DSWIS_gs(t_DSWIS_gs<0))/nd_t, d_DSWIS_gs(t_DSWIS_gs<0)/nd_u, 'LineWidth', lw);
p(3) = semilogy(abs(tD_PAR(t_PAR<0))/nd_t, d_PAR(t_PAR<0)/nd_u, 'LineWidth', lw);
p(4) = semilogy(abs(tD_EXP(t_EXP<0))/nd_t, d_EXP(t_EXP<0)/nd_u, 'LineWidth', lw);
p(5) = semilogy(abs(tD_RS(t_RS<0))/nd_t, d_RS(t_RS<0)/nd_u, 'LineWidth', lw);
p(6) = semilogy(abs(tD_SW(t_SW<0))/nd_t, d_SW(t_SW<0)/nd_u, 'LineWidth', lw);

ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$ to instability', 'Interpreter', 'latex', 'FontSize', fs);
legend([p(5) p(6) p(4) p(1) p(2) p(3) p(7)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northeast');
grid on;
ylim([1e-6 1.3]);
xlim([0 400]);

%% Energy - nondim
fig7 = figure(107);
set(gca,'FontSize',fs-2);
plot([0 0], [0 0], 'k--', 'LineWidth', lw/2);
hold on;
plot([0 0], [0 0], 'k--', 'LineWidth', lw/2);
plot([0 0], [0 0], 'k--', 'LineWidth', lw/2);
p(7) = plot([0 0 max(DSWIS_sg(:,3))]/nd_u, [1 0 0],  'k', 'LineWidth', lw);
p(1) = plot(DSWIS_sg(:,3)/nd_u, (DSWIS_sg(:,4)-mu_d*Sn)/Sn/(mu_s-mu_d), 'LineWidth', lw);
p(2) = plot(DSWIS_gs(:,3)/nd_u, (DSWIS_gs(:,4)-mu_d*Sn)/Sn/(mu_s-mu_d), 'LineWidth', lw);
p(3) = plot(PAR(:,3)/nd_u, (PAR(:,4)-mu_d*Sn)/Sn/(mu_s-mu_d), 'LineWidth', lw);
p(4) = plot(EXP(:,3)/nd_u, (EXP(:,4)-mu_d*Sn)/Sn/(mu_s-mu_d), 'LineWidth', lw);
p(5) = plot(RS(:,3)/nd_u, (RS(:,4)-mu_d*Sn)/Sn/(mu_s-mu_d), 'LineWidth', lw);
p(6) = plot(SW(:,3)/nd_u, (SW(:,4)-mu_d*Sn)/Sn/(mu_s-mu_d), 'LineWidth', lw);
p(8) = plot(EXP(:,3)/nd_u, (EXP(:,5)-mu_d*Sn)/Sn/(mu_s-mu_d), 'k--', 'LineWidth', lw);

ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend([p(5) p(6) p(4) p(1) p(2) p(3) p(7) p(8)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference', 'Loading'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northeast');
grid on;
xlim([0 1]);
ylim([-0.02 1.02]);

%% Spectra
%% Cut analytical solutions for RS
v_end     = 4*max(v_sd);

% Shear rate = 0
time_a    = time_a(Sl_r<=v_end);
Sl        = Sl(Sl_r<=v_end);
Sl_r      = Sl_r(Sl_r<=v_end);

% Shear rate = const ~= 0
time_a1   = time_a1(Sl_r1<=v_end);
Sl1       = Sl1(Sl_r1<=v_end);
Sl_r1     = Sl_r1(Sl_r1<=v_end);

vel_spec_a                              = nan(N, 9);
vel_spec_a (1:length(DSWIS_sg(:,2)),1)  = DSWIS_sg(:,2);
vel_spec_a (1:length(DSWIS_gs(:,2)),2)  = DSWIS_gs(:,2);
vel_spec_a (1:length(PAR(:,2)),3)       = PAR(:,2);
vel_spec_a (1:length(EXP(:,2)),4)       = EXP(:,2);
vel_spec_a (1:length(Sl_r1),5)          = Sl_r1;
vel_spec_a (1:length(Sl_r),6)           = Sl_r;
vel_spec_a (1:length(RS(:,2)),7)        = RS(:,2);
vel_spec_a (1:length(SW(:,2)),8)        = SW(:,2);
vel_spec_a (1:length(v_sd),9)           = v_sd;

time_spec_a                             = nan(N, 9);
time_spec_a (1:length(DSWIS_sg(:,1)),1) = DSWIS_sg(:,1);
time_spec_a (1:length(DSWIS_gs(:,1)),2) = DSWIS_gs(:,1);
time_spec_a (1:length(PAR(:,1)),3)      = PAR(:,1);
time_spec_a (1:length(EXP(:,1)),4)      = EXP(:,1);
time_spec_a (1:length(time_a1),5)       = time_a1;
time_spec_a (1:length(time_a),6)        = time_a;
time_spec_a (1:length(RS(:,1)),7)       = RS(:,1);
time_spec_a (1:length(SW(:,1)),8)       = SW(:,1);
time_spec_a (1:length(t_sd),9)          = t_sd;

sl_spec_a                               = nan(N, 9);
sl_spec_a (1:length(DSWIS_sg(:,3)),1)   = DSWIS_sg(:,3);
sl_spec_a (1:length(DSWIS_gs(:,3)),2)   = DSWIS_gs(:,3);
sl_spec_a (1:length(PAR(:,3)),3)        = PAR(:,3);
sl_spec_a (1:length(EXP(:,3)),4)        = EXP(:,3);
sl_spec_a (1:length(Sl1),5)             = Sl1;
sl_spec_a (1:length(Sl),6)              = Sl;
sl_spec_a (1:length(RS(:,3)),7)         = RS(:,3);
sl_spec_a (1:length(SW(:,3)),8)         = SW(:,3);
sl_spec_a (1:length(u_sd),9)            = u_sd;

t_sp_max    = 50*max(max(time_spec_a));
dta         = diff(time_spec_a,1);
dta(dta<=0) = NaN;
time_sp     = linspace(0,t_sp_max,200*N);
Fs          = 1/(time_sp(10)-time_sp(9));

vel_sp_pad  = zeros(length(time_sp),9);
L           = length(time_sp);
spec        = zeros(L/2,9);

% Calculate spectra
for i = 1:9
    vel_sp_pad(:,i) = interp1(time_spec_a(~isnan(time_spec_a(:,i)),i)+t_sp_max/2, vel_spec_a(~isnan(vel_spec_a(:,i)),i), time_sp, 'linear',0);
    Y               = fft(vel_sp_pad(:,i),L);
    P2              = abs(Y/L);
    P1              = P2(1:L/2+1);
    P1(2:end-1)     = 2*P1(2:end-1);
    spec(:,i)       = P1(1:L/2);
end

% Plot
fig8 = figure(108);
set(gca,'FontSize',fs-2);
loglog([0 0], [0 0], 'k--', 'LineWidth', lw/2);
hold on;
loglog([0 0], [0 0], 'k--', 'LineWidth', lw/2);
loglog([0 0], [0 0], 'k--', 'LineWidth', lw/2);
p(7) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,9)/max(spec(:,9)),  'k', 'LineWidth', lw);
p(1) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,1)/max(spec(:,9)), 'LineWidth', lw);
p(2) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,2)/max(spec(:,9)), 'LineWidth', lw);
p(3) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,3)/max(spec(:,9)), 'LineWidth', lw);
p(4) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,4)/max(spec(:,9)), 'LineWidth', lw);
p(5) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,7)/max(spec(:,9)), 'LineWidth', lw);
p(6) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,8)/max(spec(:,9)), 'LineWidth', lw);

p(10) = loglog(logspace(-1,5,100),1./(logspace(-2,4,100)), 'k:', 'LineWidth', lw/2);
p(11) = loglog(logspace(-1,5,100),1./(logspace(-2,4,100)).^2, 'k--', 'LineWidth', lw/2);
legend([p(5) p(6) p(4) p(1) p(2) p(3) p(7) p(10) p(11)],{'RS','SW','ECZ','DSWIS 1','DSWIS 2','PCZ','Reference','$f_{nd}^{\ \ -1}$','$f_{nd}^{\ \ -2}$'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$f_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlim([4e-1 0.5e2]);
ylim([5e-5 2]);
grid on;



