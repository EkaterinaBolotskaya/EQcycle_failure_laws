%% 1D spring-slider model with parabolic friction
% 03/22/2022

% by Ekaterina Bolotskaya

% This script numerically solves the equation of motion for 1D spring-slider:

% a = -1/M*(K*(D-V_0*t) + mu*sigma_n)

% where a - acceleration (second derivative of slip (D))
% M       - mass
% K       - spring stiffness
% D       - slip
% V_0     - load point velocity
% t       - time
% mu      - friction coefficient
% sigma_n - normal stress

% mu is determined by the failure law. In this case, parabolic failure law (PCZ):

% mu =
%      mu_s - (mu_s-mu_d)*((D-D_1)/D_2)^2         if D <= D_1+D_2
%      mu_d                                       if D > D_1+D_2

% where mu_s and mu_d are the static and dynamic friction coefficients, 
% D_1 and D_2 are the slip-strengthening and slip-weakening distances respectively.

clearvars;
close all;

% save data to file? 1 - yes
sd = 0;

%% Constitutive paramaters 
Sn      = 20e6;                  % normal stress at depth of interest
M       = 60e6;                  % mass

D_c     = 0.9796528896;          % characteristic slip distance, usual D_c

%% Initial and reference values
V_0     = 1e-3;                  % velocity of the load point
V_in    = 1*V_0;                 % initial sliding velocity - V_star is V*, Liu & Rice use 10-6
V_pc    = 5*V_0;                 % preseismic to coseismic phase transition
global V_ei
V_ei    = 1e-7;                  % early interseismic slip rate

%% Failure law parameters
mu_s    = 0.7;
mu_d    = 0.6; 

s_sh    = 0.117513;
s_str   = 0.735011;             % length parameters for par coh zone
mu_i    = mu_s - (mu_s-mu_d).*(-s_sh/s_str).^2;

tau_d   = mu_d*Sn;
tau_th  = mu_i*Sn;

%% Spring stiffness
K_s     = 1e6;                              % spring stiffness - factor < 1 is unstable

tf      = 2513.3*pi*sqrt(M/K_s);            % simulation length 
tpl     = tf/6;                             % plotting time

% initial slip 
u_0     = -1*mu_s*Sn/K_s;

%% Plotting and supplementary variables
lw      = 0.75;                             % line width
fs      = 11;                               % font size
mpp     = 0.05*Sn;                          % min peak prominence for peak detection

% Colors
b_col   = [0.231372549019608 0.298039215686275 0.752941176470588];
r_col   = [0.705882352941177 0.015686274509804 0.149019607843137];
g_col   = [0 .7 0];

%% Nondimentional variables for plotting
nd_t    = sqrt(M/K_s);
nd_u    = Sn*(mu_s-mu_d)/K_s;
nd_v    = nd_u/nd_t;
nd_a    = nd_v/nd_t;

%% Numerically solve 1D spring-slider equations with PAR failure law
global sl_ref_p
sl_ref_p = 0;
global k_p
k_p = 0;

% Jacobian
global sl_ref_p_j
sl_ref_p_j = 0;
global k_p_j
k_p_j = 0;

%% Non-dimentionalize
xi      = K_s*s_str/Sn;
gamma   = sqrt(K_s/M)*s_str/V_0;
Tf      = tf/(s_str/V_0);

%% dy = [du/dT dv/dT]'
dy      = @(t,y) [                        y(2)-1; 
                  -gamma*gamma*(y(1) + (1/xi)*tau_par(t, y(1), y(2), K_s, Sn, mu_s, mu_d, s_sh, s_str, V_0))];

% [u_0 v_0] initial guess
y_in    = [u_0/s_str V_0/V_0];              

% Jacobian
jac     = @(t,y) [0 1;
                 -gamma*gamma*(1+(1/xi)*jac_par(t, y(1), y(2), K_s, Sn, mu_s, mu_d, s_sh, s_str, V_0)) 0];

%% Solving a system of ordinary differential equations
options = odeset('Jacobian', jac, 'RelTol', 1e-6, 'AbsTol', 1e-6);
[T,Y]   = ode15s(dy, [0 Tf], y_in, options);

%% Original variables
v       = V_0*(Y(:,2));
t       = T*(s_str/V_0);
u       = s_str*Y(:,1);
u_abs   = V_0*t + u;
tau     = zeros(length(Y(:,1)),1);

for i = 1:length(Y(:,1))
    tau(i) = Sn*tau_par(T(i), Y(i,1), Y(i,2), K_s, Sn, mu_s, mu_d, s_sh, s_str, V_0);
end
ac      = -1/M*(K_s*u + tau);

%% Crop the time series
% Look for peaks and crossings for energy
[~,locs] = findpeaks(tau, linspace(1,length(tau),length(tau)),'MinPeakProminence',mpp);

pns     = 2;     % start from peak
time_f  = t;
tau_f   = tau;
u_f     = u;
ac_f    = ac;
v_f     = v;
u_abs_f = u_abs;
T_f     = T;

t       = t(locs(pns):end) - t(locs(pns));
t       = t(t <= tpl);
t       = [t; tpl];

tau     = tau(locs(pns):length(t)+locs(pns)-2);
tau     = [tau; interp1(time_f, tau_f, tpl+time_f(locs(pns)))];

u       = u(locs(pns):length(t)+locs(pns)-2)- u(locs(pns));
u       = [u; interp1(time_f, u_f, tpl+time_f(locs(pns)))- u_f(locs(pns))];

ac      = ac(locs(pns):length(t)+locs(pns)-2);
ac      = [ac; interp1(time_f, ac_f, tpl+time_f(locs(pns)))];

v       = v(locs(pns):length(t)+locs(pns)-2);
v       = [v; interp1(time_f, v_f, tpl+time_f(locs(pns)))];

u_abs   = u_abs(locs(pns):length(t)+locs(pns)-2)- u_abs(locs(pns));
u_abs   = [u_abs; interp1(time_f, u_abs_f, tpl+time_f(locs(pns)))- u_abs_f(locs(pns))];

T       = T(locs(pns):length(t)+locs(pns)-2);
T       = [T; interp1(time_f, T_f, tpl+time_f(locs(pns)))];

%% Vel and slip in one plot (nondim)
fig1 = figure();
set(fig1, 'defaultAxesColorOrder',[r_col; b_col]);
set(gca, 'FontSize', fs-2)
yyaxis right;
p(1) = plot(t/nd_t, u_abs/nd_u , 'LineWidth', lw);
hold on;
p(2) = plot(t/nd_t, (t*V_0+max(u))/nd_u, 'LineWidth', lw);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);

yyaxis left;
plot(t/nd_t, v/nd_v, 'LineWidth', lw);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend([p(1) p(2)],'Slider','Load point', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southeast');
set(gca,'FontSize',fs);
axis tight;
grid on;
hold off;

%% Stress Slip plot
% Look for peaks and crossings for energy
[~,locs] = findpeaks(tau+M*ac, linspace(1,length(tau),length(tau)),'MinPeakProminence',mpp);
[~,locs_t] = findpeaks(tau, linspace(1,length(tau),length(tau)),'MinPeakProminence',mpp);
lmin = find(islocalmin(tau,'MinProminence',mpp));

% Choose which peak number
pn          = 2;
id_eq       = tau <= tau_d;
u_eq        = NaN(size(u_abs));
u_eq(id_eq) = u_abs(id_eq);

[~, id_w]   = min(abs(u_eq - u_abs(locs(pn))));
u_inw       = interp1([tau(id_w) tau(id_w+1)], [u_abs(id_w) u_abs(id_w+1)], tau_d, 'linear');

id_s        = find((u_eq - u_abs(locs(pn)))>0, 1, 'first');
id_in       = locs(pn) - find(tau(locs(pn):-1:id_w) <= tau_th, 1, 'first');

idk         = locs_t(pn)+find(tau(locs_t(pn):end)+M*ac(locs_t(pn):end)<=tau_d, 1, 'first')-1;
u_ink       = interp1(tau(idk-1:idk+1)+M*ac(idk-1:idk+1), u_abs(idk-1:idk+1), tau_d);

%% Stress Slip plot (nondim)
figure()
set(gca,'FontSize', fs-2);
xb      = u_abs(locs_t(pn):id_s)/nd_u;
y1b     = (tau(locs_t(pn):id_s)-tau_d)/Sn/(mu_s-mu_d);
y2b     = (tau_d*ones(size(y1b))-tau_d)/Sn/(mu_s-mu_d);
energ_w = trapz(u_abs(locs_t(pn):id_s), tau(locs_t(pn):id_s)-tau_d);
fill([xb; flipud(xb)], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);

hold on;

xr      = [u_inw; u_abs(id_w+1:locs_t(pn))]/nd_u;
y1r     = [tau_d-tau_d; tau(id_w+1:locs_t(pn))-tau_d]/Sn/(mu_s-mu_d);
y2r     = (tau_d*ones(size(y1r))-tau_d)/Sn/(mu_s-mu_d);
energ_s = trapz([u_inw; u_abs(id_w+1:locs_t(pn))], [tau_d; tau(id_w+1:locs_t(pn))]-tau_d);
fill([xr; flipud(xr)], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);

xg      = [u_abs(locs_t(pn):idk-1); u_ink]/nd_u;
y1g     = ([tau(locs_t(pn):idk-1)+M*ac(locs_t(pn):idk-1); tau_d]-mu_d*Sn)/Sn/(mu_s-mu_d);
y2g     = (tau(locs_t(pn):idk)-mu_d*Sn)/Sn/(mu_s-mu_d);
fill([xg; flipud(xg)], [y2g; flipud(y1g)], 'k', 'EdgeColor', [.4 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.2);

p(1) = plot(u_abs/nd_u, (tau-tau_d)/Sn/(mu_s-mu_d), 'k', 'LineWidth', lw);
p(2) = plot(u_abs/nd_u, (ac*M+tau-tau_d)/Sn/(mu_s-mu_d), 'k--', 'LineWidth', lw);
legend([p(1) p(2)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
grid on;
axis tight;
hold off;
box on;

%% Phase trajectories nondimentional stress vs. ln(nondimentional velosity)
figure();
set(gca,'FontSize', fs-2);
id_pc = locs(pn)-1+find(v(locs(pn):lmin(pn+1))>=V_pc, 1, 'first');
semilogx(v(lmin(pn):locs(pn))/nd_v, (tau(lmin(pn):locs(pn)) - tau_d)/Sn/(mu_s-mu_d), 'k', 'LineWidth', lw);
hold on;
semilogx(v(id_in:locs(pn))/nd_v, (tau(id_in:locs(pn)) - tau_d)/Sn/(mu_s-mu_d), 'Color', r_col, 'LineWidth', lw);
semilogx(v(locs(pn):id_pc)/nd_v, (tau(locs(pn):id_pc) - tau_d)/Sn/(mu_s-mu_d), 'Color', g_col, 'LineWidth', lw);
semilogx(v(id_pc:lmin(pn+1))/nd_v, (tau(id_pc:lmin(pn+1)) - tau_d)/Sn/(mu_s-mu_d), 'Color', b_col, 'LineWidth', lw);

semilogx(v(lmin(pn):id_in)/nd_v, (ac(lmin(pn):id_in)*M+tau(lmin(pn):id_in)-tau_d)/Sn/(mu_s-mu_d), 'k--', 'LineWidth', lw);
semilogx(v(id_in:locs(pn))/nd_v, (ac(id_in:locs(pn))*M+tau(id_in:locs(pn))-tau_d)/Sn/(mu_s-mu_d), '--', 'Color', r_col, 'LineWidth', lw);
semilogx(v(locs(pn):id_pc)/nd_v, (ac(locs(pn):id_pc)*M+tau(locs(pn):id_pc)-tau_d)/Sn/(mu_s-mu_d), '--', 'Color', g_col, 'LineWidth', lw);
semilogx(v(id_pc:lmin(pn+1))/nd_v, (ac(id_pc:lmin(pn+1))*M+tau(id_pc:lmin(pn+1))-tau_d)/Sn/(mu_s-mu_d), '--', 'Color', b_col, 'LineWidth', lw);

semilogx(V_ei/nd_v*ones(1,100), 1.02*(linspace(min(tau(id_pc:lmin(pn+1)) - tau_d)/Sn/(mu_s-mu_d), max(ac(id_in:locs(pn))*M+tau(id_in:locs(pn))-tau_d)/Sn/(mu_s-mu_d), 100)), 'k-.', 'LineWidth', lw/2);
semilogx(V_pc/nd_v*ones(1,100), 1.02*(linspace(min(tau(id_pc:lmin(pn+1)) - tau_d)/Sn/(mu_s-mu_d), max(ac(id_in:locs(pn))*M+tau(id_in:locs(pn))-tau_d)/Sn/(mu_s-mu_d), 100)), 'k-.', 'LineWidth', lw/2);
semilogx(V_0/nd_v*ones(1,100), 1.02*(linspace(min(tau(id_pc:lmin(pn+1)) - tau_d)/Sn/(mu_s-mu_d), max(ac(id_in:locs(pn))*M+tau(id_in:locs(pn))-tau_d)/Sn/(mu_s-mu_d), 100)), 'k-.', 'LineWidth', lw/2);

text(V_0*1.1/nd_v,0.07, '$V_0$', 'Interpreter', 'latex', 'FontSize', fs-3);
text(V_pc*1.1/nd_v,0.07, '$V_{p-c}$', 'Interpreter', 'latex', 'FontSize', fs-3);
text(V_ei*1.1/nd_v,0.07, '$V_{ei}$', 'Interpreter', 'latex', 'FontSize', fs-3);

ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend({'Early interseismic', 'Late interseismic', 'Preseismic', 'Coseismic'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
grid on;
axis tight;
hold off;

%% Print
fprintf('Energies: E_s = %.4f Pa*m; ', energ_s);
fprintf('E_w = %.4f Pa*m \n', energ_w);

rec_int = t(locs(pn)) - t(locs(pn-1));
tot_sd  = tau(locs(pn)) - tau(lmin(pn-1));
fprintf('Recurrence interval = %.4f s; ', rec_int);
fprintf('Static stress drop = %.4f Pa \n', tot_sd);

fprintf('Interseismic: Slip = %.4f mm; ', (u_abs(id_in)-u_abs(lmin(pn)))*1e3);
fprintf('Velocity = %.2e ... %.2e m/s; ', min(v(lmin(pn):id_in)), max(v(lmin(pn):id_in)));
fprintf('Time = %.4f s \n', t(id_in)-t(lmin(pn)));

fprintf('Precursor: Slip = %.4f mm ', (u_abs(locs(pn))-u_abs(id_in))*1e3);
fprintf('Velocity = %.2e ... %.4f m/s ', min(v(id_in:locs(pn))), max(v(id_in:locs(pn))));
fprintf('Time = %.4f s \n', t(locs(pn))-t(id_in));

fprintf('Instability: Slip = %.4f mm; ', (u_abs(lmin(pn+1))-u_abs(locs(pn)))*1e3);
fprintf('Velocity = %.2e ... %.4f m/s; ', min(v(locs(pn):lmin(pn+1))), max(v(locs(pn):lmin(pn+1))));
fprintf('Time = %.4f s \n', t(lmin(pn+1))-t(locs(pn)));

%% Save to file
if sd == 1
    PAR = [t(lmin(pn):lmin(pn+1))-t(lmin(pn)) v(lmin(pn):lmin(pn+1)) u_abs(lmin(pn):lmin(pn+1))-u_abs(lmin(pn)) tau(lmin(pn):lmin(pn+1)) ac(lmin(pn):lmin(pn+1))*M+tau(lmin(pn):lmin(pn+1))];
    save('Data_FL/PAR.mat', 'PAR');
end

%% Function for tau parabolic
function ta = tau_par(t, u1, u2, K_s, Sn, mu_s, mu_d, slip_shift, slip_stretch, V_star) 
global sl_ref_p
global k_p
global V_ei

vel_tol = V_ei/V_star;
    if u2 > vel_tol
        if k_p == 0
            sl_ref_p = u1 + t;
            k_p = 1;
        end
        ta = -sign(u1)*max([mu_d (mu_s - (mu_s-mu_d).*(u1 - sl_ref_p + t - slip_shift/slip_stretch).^2)]);
      
    else
        ta = -sign(u1)*(mu_s - (mu_s-mu_d).*(slip_shift/slip_stretch).^2);
       
        if abs(u1*K_s*slip_stretch/Sn) <= abs(ta) && ta == -sign(u1)*(mu_s - (mu_s-mu_d).*(slip_shift/slip_stretch).^2)
            ta = -u1*K_s*slip_stretch/Sn;
            k_p = 0;
        end
    end
    
end

%% Function for Jacobian exponential
function jc = jac_par(t, u1, u2, K_s, Sn, mu_s, mu_d, slip_shift, slip_stretch, V_star) 
global sl_ref_p_j
global k_p_j
global V_ei

vel_tol = V_ei/V_star;
    if u2 > vel_tol
        if k_p_j == 0
            sl_ref_p_j = u1 + t;
            k_p_j = 1;
        end
        ta = -sign(u1)*max([mu_d (mu_s - (mu_s-mu_d).*(u1 - sl_ref_p_j + t - slip_shift/slip_stretch).^2)]);
        jc = sign(u1)*(mu_s-mu_d)*2*(u1 - sl_ref_p_j + t - slip_shift/slip_stretch);    
        if ta == -sign(u1)*mu_d
            jc = 0;
        end
    else
        ta = -sign(u1)*(mu_s - (mu_s-mu_d).*(-slip_shift/slip_stretch).^2);
        jc = sign(u1)*(mu_s-mu_d)*2*(u1 - sl_ref_p_j + t - slip_shift/slip_stretch);    
        if abs(u1*K_s*slip_stretch/Sn) <= abs(ta) && ta == -sign(u1)*(mu_d + (mu_s-mu_d).*slip_shift/slip_stretch.*exp(1-(slip_shift/slip_stretch)))
            ta = -u1*K_s*slip_stretch/Sn;
            jc = -K_s*slip_stretch/Sn;
            k_p_j = 0;
        end
    end
    
end
