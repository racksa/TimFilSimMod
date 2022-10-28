clear
close all

theta_phi_plot = figure;
theta_crit_plot = figure;

% 30 x 12
theta = [5 6.5 6.75 7 7.125 7.25 7.5 10]*pi/20;
phi = [0.058 0.058 0.058 0.058 1.53 1.53 1.53 1.53];

figure(theta_phi_plot);
hold on;
plot(theta, phi, 'k-.', 'DisplayName', '$N_x = 30, N_y = 12$');
hold off;

figure(theta_crit_plot);
hold on;
errorbar(30/12, 0.5*(7.125+7)*pi/20, 0.5*(7.125-7)*pi/20, 'k');
hold off;

% 15 x 8
theta = [5 6.5 6.625 6.75 6.875 7 10]*pi/20;
phi = [0.064 0.064 0.064 0.064 0.067 1.5046 1.5046];

figure(theta_phi_plot);
hold on;
plot(theta, phi, 'k:', 'DisplayName', '$N_x = 15, N_y = 8$');
hold off;

figure(theta_crit_plot);
hold on;
errorbar(15/8, 0.5*(7+6.875)*pi/20, 0.5*(7-6.875)*pi/20, 'k');
hold off;

% 15 x 10
theta = [5 6 6.125 6.25 6.375 6.5 10]*pi/20;
phi = [0.061 0.061 0.061 1.512 1.511 1.51 1.51];

figure(theta_phi_plot);
hold on;
plot(theta, phi, 'k-', 'DisplayName', '$N_x = 15, N_y = 10$');
hold off;

figure(theta_crit_plot);
hold on;
errorbar(15/10, 0.5*(6.25+6.125)*pi/20, 0.5*(6.25-6.125)*pi/20, 'k');
hold off;

% 15 x 14
theta = [5 5.125 5.25 5.375 5.5 10]*pi/20;
phi = [0.0559 0.0559 1.5172 1.5181 1.518 1.518];

figure(theta_phi_plot);
hold on;
plot(theta, phi, 'k--', 'DisplayName', '$N_x = 15, N_y = 14$');
hold off;

figure(theta_crit_plot);
hold on;
errorbar(15/14, 0.5*(5.125+5.25)*pi/20, 0.5*(5.25-5.125)*pi/20, 'k');
hold off;

% Figure properties
figure(theta_phi_plot);
set(gca,'FontSize',24,'FontName','Times');
xlabel('$\theta$','Interpreter','latex');
ylabel('$\phi$','Interpreter','latex');
axis([0.25*pi 0.5*pi 0 0.5*pi]);
xticks([0 pi/8 pi/4 3*pi/8 pi/2]);
yticks([0 pi/8 pi/4 3*pi/8 pi/2]);
xticklabels({'0','\pi/8','\pi/4','3\pi/8','\pi/2'});
yticklabels({'0','\pi/8','\pi/4','3\pi/8','\pi/2'});
legend1 = legend(gca,'show');
legend1.Interpreter = 'latex';
legend1.Location = 'southeast';

figure(theta_crit_plot);
set(gca,'FontSize',24,'FontName','Times');
xlabel('$N_x/N_y$','Interpreter','latex');
ylabel('$\theta_c$','Interpreter','latex');
axis([1 2.5 pi/4 6*pi/16]);
yticks([pi/4 5*pi/16 6*pi/16]);
yticklabels({'\pi/4', '5\pi/16', '6\pi/16'});

