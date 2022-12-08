clear

% The surface is curved in the y-direction and theta is measured against
% the x-axis.

theta_p = 0:pi/20:pi/2;

hold on;

% Flat
theta_b = theta_p;
plot(theta_p, theta_b, 'ko--', 'DisplayName', '$R = \infty$');

% R = -2L
theta_b = [0 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2];
plot(theta_p, theta_b, 'k^-', 'DisplayName', '$R = -2L$');

% R = -L
theta_b = [0 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2];
plot(theta_p, theta_b, 'ks-', 'DisplayName', '$R = -L$');

% R = L
theta_b = [0 0 0 0 0 0 0 0 0 0 pi/2];
plot(theta_p, theta_b, 'kd:', 'DisplayName', '$R = L$');

% R = 2L
theta_b = [0 0 0 0 0 0 0 0 0 0 pi/2];
plot(theta_p, theta_b, 'k*:', 'DisplayName', '$R = 2L$');

% Figure properties
set(gca,'FontSize',24,'FontName','Times');
xlabel('$\theta$','Interpreter','latex');
ylabel('$\phi$','Interpreter','latex');
xticks([0 pi/8 pi/4 3*pi/8 pi/2]);
yticks([0 pi/8 pi/4 3*pi/8 pi/2]);
xticklabels({'0','\pi/8','\pi/4','3\pi/8','\pi/2'});
yticklabels({'0','\pi/8','\pi/4','3\pi/8','\pi/2'});
legend1 = legend(gca,'show');
set(legend1,'Interpreter','latex');