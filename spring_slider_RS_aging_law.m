%% 1D spring-slider model with rate-and-state friction with aging law
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

% mu is determined by the failure law. In this case, rate-and-state friction (RS) with aging law:

% mu = mu_0 + a ln(V/V_0) + b ln((V_0*theta)/L)
% (theta)' = 1 - (V*theta)/L

% where mu_0 is the reference friction coefficient, V is the block slip rate, 
% theta is the state variable, representing the sliding history, V_0 is the reference slip rate, 
% L is the characteristic length scale, and a and b are rate-and-state parameters. 

clearvars;
close all;

% save data to file? 1 - yes
sd = 0;

%% Constitutive paramaters 
Sn      = 20e6;                  % normal stress at depth of interest
M       = 60e6;                  % mass

L       = 0.09;
a       = 0.0144;
b       = 0.02312;               % characteristic slip distance, usual D_c

A       = a*Sn;                       
B       = b*Sn;                        

mu_i    = 0.666519587834406;
tau_th  = mu_i*Sn;

%% Initial and reference values
V_0     = 1e-3;                 % velocity of the load point
V_pc    = 5*V_0;                % preseismic to coseismic phase transition
V_ei    = 1e-7;
V_in    = V_0;                  % initial sliding velocity - V_star is V*, Liu & Rice use 10-6
mu_s    = 0.7;
mu_d    = 0.6; 

mu_d_rs = 0.64395;

tau_st  = mu_d_rs*Sn;           % friction stress - tau_star is tau* = mu_0*Sn
tau_d   = mu_d*Sn;

%% Spring stiffness
K_s     = 1e6;
tf      = 5637*L/V_0;           % simulation length 
tpl     = 120*L/V_0;            % plotting time

% initial slip 
u_0     = -0.9*(mu_d_rs+b)*Sn/K_s;

%% Plotting and supplementary variables
lw      = 0.75;                 % line width
fs      = 11;                   % font size
mpp     = 0.8*Sn*(b-a);         % min peak prominence for peak detection

% Colors
b_col   = [0.231372549019608 0.298039215686275 0.752941176470588];
r_col   = [0.705882352941177 0.015686274509804 0.149019607843137];
g_col   = [0 .7 0];

%% Nondimentional variables for plotting
nd_t = sqrt(M/K_s);
nd_u = Sn*(mu_s-mu_d)/K_s;
nd_v = nd_u/nd_t;
nd_a = nd_v/nd_t;

%% Numerically solve 1D spring-slider equations with RS
%% Non-dimentionalize
epsilon    = (B-A)/A;
xi         = K_s*L/A;
gamma      = sqrt(K_s/M)*L/V_0;
kappa      = A*V_0/L;
tau_star_non = tau_st/K_s/L;

Tf         = tf/(L/V_0);

%% dy = [dtheta/dT du/dT dv/dT]'
dy         = @(t,y) [                                                     1/kappa-y(3)*y(1);
                                                                                     y(3)-1; 
                     -gamma*gamma*(y(2)+ (1/xi)*(mu_d_rs/a + log(y(3))+(1+epsilon)*log(kappa*y(1))))];

% [theta_0 u_0 v_0] non-dimentional initial guess   
y_in       = [13/kappa/(V_in/V_0)^(A/B) u_0/L V_in/V_0];

% Jacobian
jac        = @(t,y) [-y(3) 0 -y(1);
                         0 0     1;
                   -gamma*gamma*(1/xi)*(1+epsilon)/y(1) -gamma*gamma -gamma*gamma*(1/xi)/y(3)];

%% Solving a system of ordinary differential equations
options    = odeset('Jacobian', jac, 'RelTol', 1e-6, 'AbsTol', 1e-6);
[T,Y]      = ode15s(dy, [0 Tf], y_in, options);

%% Original variables
v          = V_0*(Y(:,3));
theta      = A*Y(:,1);
t          = T*(L/V_0);
u          = L*Y(:,2);
u_abs      = V_0*t + u;
tau        = tau_st + A*log(v/V_0) + B*log(V_0*theta/L);
mu         = tau/Sn; 

ac         = -1/M*(K_s*u + tau);

%% Crop the time series
% Look for peaks and crossings for energy
[~,locs] = findpeaks(tau, linspace(1,length(tau),length(tau)),'MinPeakProminence',mpp);

pns     = 6;     % start from peak
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

%% Vel and slip on one plot (nondim)
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
    RS = [t(lmin(pn):lmin(pn+1))-t(lmin(pn)) v(lmin(pn):lmin(pn+1)) u_abs(lmin(pn):lmin(pn+1))-u_abs(lmin(pn)) tau(lmin(pn):lmin(pn+1)) ac(lmin(pn):lmin(pn+1))*M+tau(lmin(pn):lmin(pn+1))];
    save('Data_FL/RS.mat', 'RS');
end

%% Dieterich 1992 analytical solution 
t_RS       = t(lmin(pn):lmin(pn+1))-t(lmin(pn));
v_RS       = v(lmin(pn):lmin(pn+1));
theta_RS   = theta(lmin(pn):lmin(pn+1));
ac_RS      = ac(lmin(pn):lmin(pn+1));

id_s       = find(v_RS == min(v_RS),1);
t_RS       = t_RS(id_s:end);
v_RS       = v_RS(id_s:end);
theta_RS   = theta_RS(id_s:end);
ac_RS      = ac_RS(id_s:end);


v_end      = inf;
h          = -K_s/Sn+b/L;             % term from Dieterich (model + constitutive parameters)   
v_st       = min(v(lmin(pn):lmin(pn+1)));
tau_r1     = K_s*(V_0);                                % shear rate
N          = 1e5;

% Shear rate = const = 0
t_EQ       = a/h*(1/v_st);
time_a     = fliplr(t_EQ-(logspace(-10,log10(t_EQ),N)));
Sl         = -a/h*log(1-v_st*h*time_a/a);           % slip
Sl_r       = (1/v_st-h*time_a/a).^(-1);             % slip rate

% Cut
time_a     = time_a(Sl_r<=v_end);
Sl         = Sl(Sl_r<=v_end);
Sl_r       = Sl_r(Sl_r<=v_end);

% Shear rate = const ~= 0 
t_EQ1      = a*Sn/tau_r1*log(tau_r1/h/Sn/v_st+1);     % time to instability  
time1      = fliplr(t_EQ1-(logspace(-10,log10(t_EQ1),N)));
tau1       = tau_r1*time1;
Sl1        = -a/h*log(v_st*h*Sn/tau_r1*(1-exp(tau_r1*time1/a/Sn))+1);              % slip
Sl_r1      = ((1/v_st+h*Sn/tau_r1)*exp(-tau_r1*time1/A)-h*Sn/tau_r1).^(-1);         % slip rate

% Cut
time1      = time1(Sl_r1<=v_end);
Sl1        = Sl1(Sl_r1<=v_end);
Sl_r1      = Sl_r1(Sl_r1<=v_end);


%% Dieterich plot
im_RS        = find(v_RS  == max(v_RS), 1, 'last');
im_RSa       = find(Sl_r == max(Sl_r), 1, 'last');
im_RSa1      = find(Sl_r1 == max(Sl_r1), 1, 'last');

V_int        = V_0;              % crossing for Dieterich's plot
[~,ic_RSa1]  = min(abs(Sl_r1 - V_int));
[~,ic_RSa]   = min(abs(Sl_r - V_int));
[~,ic_RS]    = min(abs(v_RS(1:im_RS) - V_int));

V_inta       = 1e6*V_0;          % crossing for Dieterich's plot 2
[~,ic_RSa1a] = min(abs(Sl_r1 - V_inta));
[~,ic_RSaa]  = min(abs(Sl_r - V_inta));

t_c          = time_a(ic_RSa);
t_m          = time_a(im_RSa);
t_ca         = time_a(ic_RSaa);
 
time_RSa1  = time1-time1(im_RSa1)+t_ca-t_m;
time_RSa   = time_a-time_a(im_RSa);
time_RS    = t_RS-t_RS(ic_RS)+t_c-t_m;

%% Plot
figure()
set(gca, 'fontsize', fs-2);
loglog(abs(time_RS(time_RS<0))/nd_t, V_0*ones(size(abs(time_RS(time_RS<0))))/nd_v, 'k-.', 'LineWidth', lw/2);
hold on;
loglog(abs(time_RS(time_RS<0))/nd_t, V_ei*ones(size(abs(time_RS(time_RS<0))))/nd_v, 'k-.', 'LineWidth', lw/2);
loglog(abs(time_RS(time_RS<0))/nd_t, V_pc*ones(size(abs(time_RS(time_RS<0))))/nd_v, 'k-.', 'LineWidth', lw/2);

p(1)=loglog(abs(time_RSa(time_RSa<0))/nd_t, Sl_r(time_RSa<0)/nd_v, 'LineWidth', lw);
p(2)=loglog(abs(time_RSa1(time_RSa1<0))/nd_t, Sl_r1(time_RSa1<0)/nd_v, '--', 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', lw);
p(3)=loglog(abs(time_RS(time_RS<0))/nd_t, v_RS(time_RS<0)/nd_v, 'LineWidth', lw);
p(4)=loglog(abs(time_RS(time_RS<0))/nd_t, L./theta_RS(time_RS<0)/nd_v, 'LineWidth', lw);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$\ to\ instability', 'Interpreter', 'latex', 'FontSize', fs);
legend(p(1:4),{'RS analytical $\dot\tau = 0$','RS analytical $\dot\tau \neq 0$', 'RS\ spring-slider','L/$\theta$'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
text(100,V_0*1.5/nd_v, '$V_0$', 'Interpreter', 'latex', 'FontSize', fs-1);
text(100,V_pc*1.5/nd_v, '$V_{p-c}$', 'Interpreter', 'latex', 'FontSize', fs-1);
text(100,V_ei*1.5/nd_v, '$V_{ei}$', 'Interpreter', 'latex', 'FontSize', fs-1);

grid on;
axis tight;
xlim([0.042232 max(t_RS)/nd_t]);
ylim([0.5*min(v_RS/nd_v) 1]);

