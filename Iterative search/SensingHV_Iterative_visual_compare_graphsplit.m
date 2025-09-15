clc
clear
close all;

%rng(1); % Random 상수 고정

%% 시뮬레이션 상황 변수 선언

P = 4; % N of NLoS path

SV = [0; 0; 0]; % Sensing Vehicle position
HV = [8; -20; 5]; % Hidden Vehicle position
sc1 = [5; -5; 5+4*rand(1)]; % Scatterer 1 
sc2 = [15; -5; 5+4*rand(1)]; % Scatterer 2
sc3 = [15; 5; 5+4*rand(1)]; % Scatterer 3
sc4 = [5; 5; 5+4*rand(1)]; % Scatterer 4

c = 3e8; % speed of light
Q_true = deg2rad(5); % elevation of HV (answer)
w_true = deg2rad(70); % angle of HV compare to SV (answer)

% Z angle of AoA
alpha = [atan2(sqrt(sc1(1)^2+sc1(2)^2),sc1(3)); atan2(sqrt(sc2(1)^2+sc2(2)^2),sc2(3)); atan2(sqrt(sc3(1)^2+sc3(2)^2),sc3(3)); atan2(sqrt(sc4(1)^2+sc4(2)^2),sc4(3))];
% AoA
theta = [atan2(sc1(2),sc1(1)); atan2(sc2(2),sc2(1)); atan2(sc3(2),sc3(1)); atan2(sc4(2),sc4(1))];
% Z angle of AoD
psipos = [atan2(sqrt((sc1(1)-HV(1))^2+(sc1(2)-HV(2))^2),sc1(3)-HV(3)); atan2(sqrt((sc2(1)-HV(1))^2+(sc2(2)-HV(2))^2),sc2(3)-HV(3)); atan2(sqrt((sc3(1)-HV(1))^2+(sc3(2)-HV(2))^2),sc3(3)-HV(3)); atan2(sqrt((sc4(1)-HV(1))^2+(sc4(2)-HV(2))^2),sc4(3)-HV(3))];
psi = psipos - Q_true;
% AoD
phipos = [atan2(sc1(2)-HV(2),sc1(1)-HV(1)); atan2(sc2(2)-HV(2),sc2(1)-HV(1)); atan2(sc3(2)-HV(2),sc3(1)-HV(1)); atan2(sc4(2)-HV(2),sc4(1)-HV(1))];
phi = phipos - w_true;

% d1, TDoA
d1 = sqrt(sc1(1)^2+sc1(2)^2+sc1(3)^2)+sqrt((sc1(1)-HV(1))^2+(sc1(2)-HV(2))^2+(sc1(3)-HV(3))^2);
d2 = sqrt(sc2(1)^2+sc2(2)^2+sc2(3)^2)+sqrt((sc2(1)-HV(1))^2+(sc2(2)-HV(2))^2+(sc2(3)-HV(3))^2);
d3 = sqrt(sc3(1)^2+sc3(2)^2+sc3(3)^2)+sqrt((sc3(1)-HV(1))^2+(sc3(2)-HV(2))^2+(sc3(3)-HV(3))^2);
d4 = sqrt(sc4(1)^2+sc4(2)^2+sc4(3)^2)+sqrt((sc4(1)-HV(1))^2+(sc4(2)-HV(2))^2+(sc4(3)-HV(3))^2);
tdoa = [(d1-d1)/c;(d2-d1)/c;(d3-d1)/c;(d4-d1)/c];

%%% Adding noise %%%
% Random noise per iteration
% alpha = alpha + deg2rad(1) * randn(size(alpha));
% theta = theta + deg2rad(1) * randn(size(theta));
% psi = psi + deg2rad(1) * randn(size(psi));
% phi = phi + deg2rad(1) * randn(size(phi));
% tdoa = tdoa + 1e-9 * randn(size(tdoa))* 1;
% d1 = d1 + 0.1 * randn(1);

%% 비교를 위한 위치 단일 추정 알고리즘
Qpos = deg2rad(0);
wpos = deg2rad(90);

vpos = calvp(Qpos,wpos,alpha,theta,psi,phi,tdoa,d1,c);

xpos = vpos.*(sin(alpha).*cos(theta)+sin(psi).*cos(phi))-(d1+c*tdoa).*sin(psi).*cos(phi);
ypos = vpos.*(sin(alpha).*sin(theta)+sin(psi).*sin(phi))-(d1+c*tdoa).*sin(psi).*sin(phi);
zpos = vpos.*(cos(alpha)+cos(psi))-(d1+c*tdoa).*cos(psi);

x_estimpos = sum(vpos.*(sin(alpha).*cos(theta)+sin(psi+Qpos).*cos(phi+wpos))-(d1+c*tdoa).*sin(psi+Qpos).*cos(phi+wpos))/size(vpos,1)
y_estimpos = sum(vpos.*(sin(alpha).*sin(theta)+sin(psi+Qpos).*sin(phi+wpos))-(d1+c*tdoa).*sin(psi+Qpos).*sin(phi+wpos))/size(vpos,1)
z_estimpos = sum(vpos.*(cos(alpha)+cos(psi+Qpos))-(d1+c*tdoa).*cos(psi+Qpos))/size(vpos,1)

%% NLoS 기반 추정 알고리즘

Q0 = deg2rad(1); % Q 초기값
w0 = deg2rad(90); % w 초기값
v0 = calvp(Q0,w0,alpha,theta,psi,phi,tdoa,d1,c); % v 초기값 행렬

Q_prev = Q0; % 직전 iteration 값 설정 : Q(t-1)
w_prev = w0; % 직전 iteration 값 설정 : w(t-1)
v_prev = v0; % 직전 iteration 행렬 값 설정 : v(t-1)

tol = 3e-3; % 수렴 조건
deltaQ = 1; % 수렴 확인을 위한 변수 (iteration에 따른 Q 변화량)
deltaw = 1; % 수렴 확인을 위한 변수 (iteration에 따른 w 변화량)
deltav = ones(P,1); % 수렴 확인을 위한 변수 (iteration에 따른 v 변화량)

count = 0; % 수렴 loop 횟수 기록용

% Q_state와 w_state 값을 기록할 벡터 초기화
Q_values = []; 
w_values = [];
v_values = [];

iterations = [];
tic;
while (abs(rad2deg(deltaQ)) > tol) && (abs(rad2deg(deltaw)) > tol) && all(abs(deltav) > tol)

count = count+1;

% A*X=B
A = matA_caliter(psi, phi, Q_prev, w_prev, tdoa, d1, v_prev, c);
B = matB_caliter(alpha, theta, psi, phi, Q_prev, w_prev, tdoa, d1, v_prev, c);
% Calculating X (X=[Q;w])

X = A\B; % Psuedo-Inverse

deltaQ = X(1);
deltaw = X(2);

Q_state = deltaQ+Q_prev;
w_state = deltaw+w_prev;
v_state = calvp(Q_state,w_state,alpha,theta,psi,phi,tdoa,d1,c);

Q_prev = Q_state;
w_prev = w_state;
deltav = v_state-v_prev;
v_prev = v_state;

% Q_state와 w_state를 기록
Q_values = [Q_values, Q_state];
w_values = [w_values, w_state];
v_values = [v_values, v_state];
iterations = [iterations, count];

fprintf('Iteration : %d, Q(t) : %f , w(t) : %f\n',count,rad2deg(Q_state),rad2deg(w_state));

end
elapsed = toc;
avg_iter_time = elapsed / count;
fprintf('계산 소요 시간 : %.6f ms\n ,평균 1회 iteration 시간: %.6f ms\n', elapsed*1000 , avg_iter_time*1000);

x_estim0 = sum(v0.*(sin(alpha).*cos(theta)+sin(psi+Q0).*cos(phi+w0))-(d1+c*tdoa).*sin(psi+Q0).*cos(phi+w0))/size(v0,1);
y_estim0 = sum(v0.*(sin(alpha).*sin(theta)+sin(psi+Q0).*sin(phi+w0))-(d1+c*tdoa).*sin(psi+Q0).*sin(phi+w0))/size(v0,1);
z_estim0 = sum(v0.*(cos(alpha)+cos(psi+Q0))-(d1+c*tdoa).*cos(psi+Q0))/size(v0,1);

for i = 1:size(iterations,2)
x_est = v_values(:,i).*(sin(alpha).*cos(theta)+sin(psi+Q_values(i)).*cos(phi+w_values(i)))-(d1+c*tdoa).*sin(psi+Q_values(i)).*cos(phi+w_values(i));
y_est = v_values(:,i).*(sin(alpha).*sin(theta)+sin(psi+Q_values(i)).*sin(phi+w_values(i)))-(d1+c*tdoa).*sin(psi+Q_values(i)).*sin(phi+w_values(i));
z_est = v_values(:,i).*(cos(alpha)+cos(psi+Q_values(i)))-(d1+c*tdoa).*cos(psi+Q_values(i));

x_estim(i) = sum(x_est)/size(x_est,1);
y_estim(i) = sum(y_est)/size(y_est,1);
z_estim(i) = sum(z_est)/size(z_est,1);
end

% while문 종료 후, subplot을 사용하여 Q_state와 w_state의 변화 시각화
t = tiledlayout(2,2);

% 왼쪽(1,1)에 큰 플롯 하나
nexttile(1, [2 1]); % 2행 1열을 합쳐서 사용
hold on;
grid on;

walls = [
    -5, 5, 0; -5, 5, 15; 5, 5, 15; 5, 5, 0;
    15, 5, 0; 15, 5, 15; 25, 5, 15; 25, 5, 0;
    5, 20, 0; 5, 20, 15; 5, 5, 15; 5, 5, 0;
    15, 20, 0; 15, 20, 15; 15, 5, 15; 15, 5, 0;
    -5, -5, 0; -5, -5, 15; 5, -5, 15; 5, -5, 0;
    15, -5, 0; 15, -5, 15; 15, -30, 15; 15, -30, 0;
    15, -5, 0; 15, -5, 15; 25, -5, 15; 25, -5, 0;
    5, -5, 0; 5, -5, 15; 5, -30, 15; 5, -30, 0
]; % visualize walls

grounds = [
    -5, 5, 0; -5, -5, 0; 18, -5, 0; 18, 5, 0;
    18, -5, 0; 18, 5, 0; 25, 5, 3; 25, -5, 3;
    5, -30, 4.1; 5, -15, 5.3; 15, -15, 5.3; 15, -30, 4.1;
    5, -15, 5.3; 5, -5, 0; 15, -5, 0; 15, -15, 5.3;
    5, -30, 4.1; 5, -30, 0; 5, -5, 0; 5, -15, 5.3;
    15, -30, 4.1; 15, -30, 0; 15, -5, 0; 15, -15, 5.3;
    15, 20, 0; 5, 20, 0; 5, 5, 0; 15, 5, 0
]; % visualize grounds

scatter3(SV(1), SV(2), SV(3), 100, 'r', 'filled', 'DisplayName', 'SV');
scatter3(HV(1), HV(2), HV(3), 100, 'b', 'filled', 'DisplayName', 'HV');

scatter3(x_estimpos, y_estimpos, z_estimpos, 100, 'cyan', 'filled', 'DisplayName', 'HV');

scatter3(x_estim(size(iterations,2)), y_estim(size(iterations,2)), z_estim(size(iterations,2)), 150, 'm', 'filled', 'DisplayName', 'Estimation');

scatter3(sc1(1), sc1(2), sc1(3), 100, 'k', 'filled', 'DisplayName', 'Scatterer 1');
scatter3(sc2(1), sc2(2), sc2(3), 100, 'k', 'filled', 'DisplayName', 'Scatterer 2');
scatter3(sc3(1), sc3(2), sc3(3), 100, 'k', 'filled', 'DisplayName', 'Scatterer 3');
scatter3(sc4(1), sc4(2), sc4(3), 100, 'k', 'filled', 'DisplayName', 'Scatterer 4');

plot3([SV(1) sc1(1)], [SV(2) sc1(2)], [SV(3) sc1(3)], 'k--', 'LineWidth', 1.5);
plot3([sc1(1) HV(1)], [sc1(2) HV(2)], [sc1(3) HV(3)], 'k--', 'LineWidth', 1.5);

plot3([SV(1) sc2(1)], [SV(2) sc2(2)], [SV(3) sc2(3)], 'k--', 'LineWidth', 1.5);
plot3([sc2(1) HV(1)], [sc2(2) HV(2)], [sc2(3) HV(3)], 'k--', 'LineWidth', 1.5);

plot3([SV(1) sc3(1)], [SV(2) sc3(2)], [SV(3) sc3(3)], 'k--', 'LineWidth', 1.5);
plot3([sc3(1) HV(1)], [sc3(2) HV(2)], [sc3(3) HV(3)], 'k--', 'LineWidth', 1.5);

plot3([SV(1) sc4(1)], [SV(2) sc4(2)], [SV(3) sc4(3)], 'k--', 'LineWidth', 1.5);
plot3([sc4(1) HV(1)], [sc4(2) HV(2)], [sc4(3) HV(3)], 'k--', 'LineWidth', 1.5);

quiver3(SV(1), SV(2), SV(3), 7, 0, 0, 0, 'r', 'LineWidth', 2);
[dx_ans, dy_ans, dz_ans] = sph2cart(w_true, Q_true, 7);
quiver3(HV(1), HV(2), HV(3), dx_ans, dy_ans, dz_ans, 0, 'b', 'LineWidth', 2);
[dx_estim, dy_estim, dz_estim] = sph2cart(w_state, Q_state, 8);
quiver3(x_estim(size(iterations,2)), y_estim(size(iterations,2)), z_estim(size(iterations,2)), dx_estim, dy_estim, dz_estim, 0, 'm', 'LineWidth', 2.5);

fill3(walls(1:4,1), walls(1:4,2), walls(1:4,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(5:8,1), walls(5:8,2), walls(5:8,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(9:12,1), walls(9:12,2), walls(9:12,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(13:16,1), walls(13:16,2), walls(13:16,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(17:20,1), walls(17:20,2), walls(17:20,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(21:24,1), walls(21:24,2), walls(21:24,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(25:28,1), walls(25:28,2), walls(25:28,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(walls(29:32,1), walls(29:32,2), walls(29:32,3), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

fill3(grounds(1:4,1), grounds(1:4,2), grounds(1:4,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(grounds(5:8,1), grounds(5:8,2), grounds(5:8,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(grounds(9:12,1), grounds(9:12,2), grounds(9:12,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(grounds(13:16,1), grounds(13:16,2), grounds(13:16,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(grounds(17:20,1), grounds(17:20,2), grounds(17:20,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

fill3(grounds(21:24,1), grounds(21:24,2), grounds(21:24,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(grounds(25:28,1), grounds(25:28,2), grounds(25:28,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

title('Simulation of estimating information of HV','FontSize',15);
legend('SV', 'HV', 'Prev','Estimation','Scatterer','','','','NLoS Path','FontSize',12);
hold off;

% 오차 백분율 계산
Q_error_percent = abs(rad2deg([Q0, Q_values]) - rad2deg(Q_true)) / abs(rad2deg(Q_true)) * 100;
w_error_percent = abs(rad2deg([w0, w_values]) - rad2deg(w_true)) / abs(rad2deg(w_true)) * 100;

% 첫 번째 서브플롯: Q_state 변화
% 오른쪽 위(1,2)
nexttile(2);

yyaxis left; % 왼쪽 y축 활성화
plot([0, iterations], [rad2deg(Q0), rad2deg(Q_values)], 'LineWidth',2,'Marker','|', 'DisplayName', 'Q State'); % Q_state의 변화
hold on;
plot([0, iterations], repmat(rad2deg(Q_true), 1, length(iterations)+1),'r--', 'LineWidth',2, 'DisplayName', 'True Q'); % Q_true의 직선 표시
y_range_left = [min(rad2deg(Q_values)) - 10, max(rad2deg(Q_values)) + 10]; % y축 범위 설정
ylim(y_range_left);
xlabel('Iteration','FontSize',15);
xticks(0:1:max(iterations));
ylabel('Q State (degrees)','FontSize',15);
title('Estimation of HV elevation State over Iterations','FontSize',15);

yyaxis right; % 오른쪽 y축 활성화
plot([0, iterations], Q_error_percent, 'LineWidth',1.5,'Marker','o','LineStyle',':', 'DisplayName', 'Q Error (%)'); % Q 오차 백분율 플롯
ylabel('Error (%)','FontSize',15);
y_range_right_q = [0, max(Q_error_percent) + 10]; % 오른쪽 y축 범위 설정 (필요시 조정)
ylim(y_range_right_q);

legend('show','Location','best','FontSize',10); % 범례 표시 (위치 및 폰트 크기 조정)
hold off;

% 두 번째 서브플롯: w_state 변화
% 오른쪽 아래(2,2)
nexttile(4);

yyaxis left; % 왼쪽 y축 활성화
plot([0, iterations],[rad2deg(w0), rad2deg(w_values)], 'LineWidth',2,'Marker','|', 'DisplayName', 'w State'); % w_state의 변화
hold on;
plot([0, iterations], repmat(rad2deg(w_true), 1, length(iterations)+1), 'r--', 'LineWidth',2, 'DisplayName', 'True w'); % w_true의 직선 표시
y_range_left_w = [min(rad2deg(w_values)) - 10, max(rad2deg(w_values)) + 20]; % y축 범위 설정
ylim(y_range_left_w);
xlabel('Iteration','FontSize',15);
xticks(0:1:max(iterations));
ylabel('w State (degrees)','FontSize',15);
title('Estimation of HV azimuth State over Iterations','FontSize',15);

yyaxis right; % 오른쪽 y축 활성화
plot([0, iterations], w_error_percent, 'LineWidth',1.5,'Marker','o','LineStyle',':', 'DisplayName', 'w Error (%)'); % w 오차 백분율 플롯
ylabel('Error (%)','FontSize',15);
y_range_right_w = [0, max(w_error_percent) + 10]; % 오른쪽 y축 범위 설정 (필요시 조정)
ylim(y_range_right_w);

legend('show','Location','best','FontSize',10); % 범례 표시 (위치 및 폰트 크기 조정)
hold off;


% 
% figure;
% subplot(1,2,1)
% hold on;
% plot([0, iterations], [100*abs(rad2deg(Q_true)-rad2deg(Q0))/rad2deg(Q_true), 100*abs(rad2deg(Q_true)-rad2deg(Q_values))./rad2deg(Q_true)], 'LineWidth',2,'Marker','|'); % Q_state 오차율의 변화
% plot([0, iterations], [100*abs(rad2deg(w_true)-rad2deg(w0))/rad2deg(w_true), 100*abs(rad2deg(w_true)-rad2deg(w_values))./rad2deg(w_true)], 'LineWidth',2,'Marker','|'); % w_state 오차율의 변화
% xlabel('Iteration','FontSize',15);
% xticks(0:1:size(iterations,2));
% ylabel('Error (%)','FontSize',15);
% title('Error of Angle State over Iterations','FontSize',15);
% legend('Q State Error', 'w State Error','FontSize',12);
% 
% hold off;
% subplot(1,2,2)
% hold on;
% plot([0, iterations], [100*abs(HV(1)-x_estim0)/abs(HV(1)), 100*abs(HV(1)-x_estim)./abs(HV(1))], 'LineWidth',2,'Marker','|'); % x_state 오차율의 변화
% plot([0, iterations], [100*abs(HV(2)-y_estim0)/abs(HV(2)), 100*abs(HV(2)-y_estim)./abs(HV(2))], 'LineWidth',2,'Marker','|'); % y_state 오차율의 변화
% plot([0, iterations], [100*abs(HV(3)-z_estim0)/abs(HV(3)), 100*abs(HV(3)-z_estim)./abs(HV(3))], 'LineWidth',2,'Marker','|'); % z_state 오차율의 변화
% xlabel('Iteration','FontSize',15);
% xticks(0:1:size(iterations,2));
% ylabel('Error (%)','FontSize',15);
% title('Error of Position State over Iterations','FontSize',15);
% legend('x State Error', 'y State Error', 'z State Error','FontSize',12);
% hold off;
