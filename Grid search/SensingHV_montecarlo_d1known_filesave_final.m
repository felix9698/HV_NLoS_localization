clc
clear
close all;

rng(1);

%%% Variables %%%

P = 3; % N of NLoS path

alpha = zeros(P,1).'; % Z angle of AoA
theta = zeros(P,1).'; % AoA
psi = zeros(P,1).'; % Z angle of AoD
phi = zeros(P,1).'; % AoD
tdoa = zeros(P,1).'; % TDoA

c = 3e8; % speed of light
Q = 0; % elevation of HV
omega = 0; % angle of HV compare to SV

numiter = 100; % N of iterations

%%% calculating variables %%%
% using this section in simulation, and change these values in real time experiment
x_monte_avg_M = []; 
y_monte_avg_M = []; 
z_monte_avg_M = []; 

for HVx = 6:2:14
for HVy = -50:5:-5

SV = [0; 0; 0]; % Sensing Vehicle
HV = [HVx; HVy; 5]; % Hidden Vehicle
sc1 = [5; -5; 3]; % Scatterer 1 
sc2 = [15; -5; 2]; % Scatterer 2
sc3 = [15; 5; 8]; % Scatterer 3

Q = deg2rad(-3); % elevation of HV
omega = deg2rad(90); % angle of HV compare to SV
 
alpha = [atan2(sqrt(sc1(1)^2+sc1(2)^2),sc1(3)); atan2(sqrt(sc2(1)^2+sc2(2)^2),sc2(3)); atan2(sqrt(sc3(1)^2+sc3(2)^2),sc3(3))];
theta = [atan2(sc1(2),sc1(1)); atan2(sc2(2),sc2(1)); atan2(sc3(2),sc3(1))];
psi = [atan2(sqrt((sc1(1)-HV(1))^2+(sc1(2)-HV(2))^2),sc1(3)-HV(3)); atan2(sqrt((sc2(1)-HV(1))^2+(sc2(2)-HV(2))^2),sc2(3)-HV(3)); atan2(sqrt((sc3(1)-HV(1))^2+(sc3(2)-HV(2))^2),sc3(3)-HV(3))] - Q;
phi = [atan2(sc1(2)-HV(2),sc1(1)-HV(1)); atan2(sc2(2)-HV(2),sc2(1)-HV(1)); atan2(sc3(2)-HV(2),sc3(1)-HV(1))] - omega;

d1 = sqrt(sc1(1)^2+sc1(2)^2+sc1(3)^2)+sqrt((sc1(1)-HV(1))^2+(sc1(2)-HV(2))^2+(sc1(3)-HV(3))^2);
d2 = sqrt(sc2(1)^2+sc2(2)^2+sc2(3)^2)+sqrt((sc2(1)-HV(1))^2+(sc2(2)-HV(2))^2+(sc2(3)-HV(3))^2);
d3 = sqrt(sc3(1)^2+sc3(2)^2+sc3(3)^2)+sqrt((sc3(1)-HV(1))^2+(sc3(2)-HV(2))^2+(sc3(3)-HV(3))^2);
tdoa = [(d1-d1)/c;(d2-d1)/c;(d3-d1)/c];

x_monte = zeros(1,numiter);
y_monte = zeros(1,numiter);
z_monte = zeros(1,numiter);

omega_monte = zeros(1,numiter);
Q_monte = zeros(1,numiter);

parfor monte = 1:numiter

%%% Adding noise %%%
% Random noise per iteration
alpha_err = alpha + deg2rad(1) * randn(size(alpha));
theta_err = theta + deg2rad(1) * randn(size(theta));
psi_err = psi + deg2rad(1) * randn(size(psi));
phi_err = phi + deg2rad(1) * randn(size(phi));
tdoa_err = tdoa + 1e-9 * randn(size(tdoa))* 1;
d1_err = d1 + 0.1 * randn(1);
 
%%% Localization algorithm %%%
w_est = [];

for omega_iter = -179:1:180
    for Q_iter = -9:1:10
        %%% matrix A %%%
        A = matA_d1known(alpha_err, theta_err, psi_err, phi_err, deg2rad(omega_iter), deg2rad(Q_iter));
        %%% matrix B %%%
        B = matB_d1known(c, tdoa_err, psi_err, deg2rad(Q_iter), deg2rad(omega_iter), phi_err, d1);

        w_est(omega_iter+180,Q_iter+10) = sum(null(A.').'*B);
    end
end
% A(w)*z = B(w)

[~,ind] = min(abs(w_est),[],'all');
Q_est = deg2rad(ceil(ind/360)-10);
omega_est = deg2rad(mod(ind,360)-180);

%fprintf('W = %.0f°, Q = %.0f°\n', rad2deg(omega_est), rad2deg(Q_est));

A = matA_d1known(alpha_err, theta_err, psi_err, phi_err, omega_est, Q_est);
B = matB_d1known(c, tdoa_err, psi_err, Q_est, omega_est, phi_err, d1);

% z = (v1,...,vp)'
z_est = inv(A.'*A)*A.'*B;
x_est_arr = z_est(1:P).*(sin(alpha_err).*cos(theta_err)+sin(psi_err+Q_est).*cos(phi_err+omega_est))-(d1+c*tdoa_err).*sin(psi_err+Q_est).*cos(phi_err+omega_est);
y_est_arr = z_est(1:P).*(sin(alpha_err).*sin(theta_err)+sin(psi_err+Q_est).*sin(phi_err+omega_est))-(d1+c*tdoa_err).*sin(psi_err+Q_est).*sin(phi_err+omega_est);
z_est_arr = z_est(1:P).*(cos(alpha_err)+cos(psi_err+Q_est))-(d1+c*tdoa_err).*cos(psi_err+Q_est);

x_estim = sum(x_est_arr)/size(x_est_arr,1);
y_estim = sum(y_est_arr)/size(y_est_arr,1);
z_estim = sum(z_est_arr)/size(z_est_arr,1);

x_monte(1,monte) = x_estim; 
y_monte(1,monte) = y_estim; 
z_monte(1,monte) = z_estim; 

omega_monte(1,monte) = omega_est; 
Q_monte(1,monte) = Q_est; 

end

x_monte_avg = rmse(x_monte,HV(1));
y_monte_avg = rmse(y_monte,HV(2));
z_monte_avg = rmse(z_monte,HV(3));

x_monte_avg_M = [x_monte_avg_M; [HVx, HVy, x_monte_avg]]; 
y_monte_avg_M = [y_monte_avg_M; [HVx, HVy, y_monte_avg]]; 
z_monte_avg_M = [z_monte_avg_M; [HVx, HVy, z_monte_avg]]; 

%%
%%% Visualization %%%
resfig = figure;
hold on;
scatter3(x_monte, y_monte, z_monte)

walls = [
    -5, 5, 0; -5, 5, 10; 25, 5, 10; 25, 5, 0;   
    -5, -5, 0; -5, -5, 10; 5, -5, 10; 5, -5, 0;
    15, -5, 0; 15, -5, 10; 25, -5, 10; 25, -5, 0;
    5, -5, 0; 5, -5, 10; 5, -50, 10; 5, -50, 0;
    15, -5, 0; 15, -5, 10; 15, -50, 10; 15, -50, 0
]; % visualize walls

%figure;
hold on;
grid on;

scatter3(SV(1), SV(2), SV(3), 100, 'r', 'filled', 'DisplayName', 'SV');
scatter3(HV(1), HV(2), HV(3), 100, 'b', 'filled', 'DisplayName', 'HV');

scatter3(sc1(1), sc1(2), sc1(3), 100, 'g', 'DisplayName', 'Scatterer 1');
scatter3(sc2(1), sc2(2), sc2(3), 100, 'g', 'DisplayName', 'Scatterer 2');
scatter3(sc3(1), sc3(2), sc3(3), 100, 'g', 'DisplayName', 'Scatterer 3');

plot3([SV(1) sc1(1)], [SV(2) sc1(2)], [SV(3) sc1(3)], 'r-', 'LineWidth', 1.5);
plot3([sc1(1) HV(1)], [sc1(2) HV(2)], [sc1(3) HV(3)], 'r--', 'LineWidth', 1.5);

plot3([SV(1) sc2(1)], [SV(2) sc2(2)], [SV(3) sc2(3)], 'r-', 'LineWidth', 1.5);
plot3([sc2(1) HV(1)], [sc2(2) HV(2)], [sc2(3) HV(3)], 'r--', 'LineWidth', 1.5);

plot3([SV(1) sc3(1)], [SV(2) sc3(2)], [SV(3) sc3(3)], 'r', 'LineWidth', 1.5);
plot3([sc3(1) HV(1)], [sc3(2) HV(2)], [sc3(3) HV(3)], 'r--', 'LineWidth', 1.5);

fill3(walls(1:4,1), walls(1:4,2), walls(1:4,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(5:8,1), walls(5:8,2), walls(5:8,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(9:12,1), walls(9:12,2), walls(9:12,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(13:16,1), walls(13:16,2), walls(13:16,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(17:20,1), walls(17:20,2), walls(17:20,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('X');
ylabel('Y');
zlabel('Z');

quiver3(SV(1), SV(2), SV(3), 6, 0, 0, 0, 'k', 'LineWidth', 2);
[dx_ans, dy_ans, dz_ans] = sph2cart(omega, Q, 8);
quiver3(HV(1), HV(2), HV(3), dx_ans, dy_ans, dz_ans, 0, 'k', 'LineWidth', 2);

[dx, dy, dz] = sph2cart(omega_monte, Q_monte, 3);
quiver3(x_monte, y_monte, z_monte, dx, dy, dz, 0, 'm', 'LineWidth', 2);

hold off;

fprintf('Now saving... HVx = %d, HVy = %d\n', HVx, HVy);

resultFolder = './result';
fileName = sprintf('err1_HVx=%d_HVy=%d.fig', HVx,HVy);
fullPath = fullfile(resultFolder, fileName); 
savefig(resfig, fullPath); 

close(resfig);

errhist = xyzerr(HV.',[x_monte.',y_monte.',z_monte.']);

fileName = sprintf('err1_HVx=%d_HVy=%d_hist.fig', HVx,HVy);
fullPath = fullfile(resultFolder, fileName); 
savefig(errhist, fullPath); 

close(errhist);

end
end
%%
figure;

hold on;
grid on;

scatter3(sc1(1), sc1(2), sc1(3), 100, 'g', 'DisplayName', 'Scatterer 1');
scatter3(sc2(1), sc2(2), sc2(3), 100, 'g', 'DisplayName', 'Scatterer 2');
scatter3(sc3(1), sc3(2), sc3(3), 100, 'g', 'DisplayName', 'Scatterer 3');

fill3(walls(1:4,1), walls(1:4,2), walls(1:4,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(5:8,1), walls(5:8,2), walls(5:8,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(9:12,1), walls(9:12,2), walls(9:12,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(13:16,1), walls(13:16,2), walls(13:16,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(17:20,1), walls(17:20,2), walls(17:20,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('X');
ylabel('Y');
zlabel('Z');

sx = scatter3(x_monte_avg_M(:,1) ,x_monte_avg_M(:,2), x_monte_avg_M(:,3),'r','filled');
sy = scatter3(y_monte_avg_M(:,1) ,y_monte_avg_M(:,2), y_monte_avg_M(:,3),'y','filled');
sz = scatter3(z_monte_avg_M(:,1) ,z_monte_avg_M(:,2), z_monte_avg_M(:,3),'b','filled');

legend([sx sy sz],{'err avg of x','err avg of y','err avg of z'})

hold off;

%%
x = x_monte_avg_M(:,1);  % x좌표 (6부터 14까지 2 단위)
y = x_monte_avg_M(:,2);  % y좌표 (-50부터 -5까지 5 단위)
values = (x_monte_avg_M(:,3)+y_monte_avg_M(:,3)+z_monte_avg_M(:,3))/3;  % 각 격자점의 값

% x와 y의 고유한 값을 얻습니다.
x_unique = unique(x);
y_unique = unique(y);

% x와 y에 대한 그리드를 생성합니다.
[X_grid, Y_grid] = meshgrid(x_unique, y_unique);

% values를 그리드에 맞게 행렬로 재구성합니다.
% Z는 y축이 행(row), x축이 열(column)로 구성된 행렬입니다.
Z = zeros(length(y_unique), length(x_unique));

% 각 데이터 포인트에 대해 Z 행렬에 값을 채웁니다.
for k = 1:length(values)
    xi = find(x_unique == x(k));
    yi = find(y_unique == y(k));
    Z(yi, xi) = values(k);
end

% 컬러맵을 그립니다.
figure;
imagesc(x_unique, y_unique, Z);
set(gca, 'YDir', 'normal');  % y축이 아래에서 위로 증가하도록 설정
colorbar;
xlabel('X coordinate of HV');
ylabel('Y coordinate of HV');
title('RMSE average');

colormap(flipud(parula));