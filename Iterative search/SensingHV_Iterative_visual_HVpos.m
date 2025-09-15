clc
clear
close all;

%rng(1); % Random 상수 고정

%% 시뮬레이션 상황 변수 선언

P = 4; % N of NLoS path

sc1 = [5; -5; 5+4*rand(1)]; % Scatterer 1 
sc2 = [15; -5; 5+4*rand(1)]; % Scatterer 2
sc3 = [15; 5; 5+4*rand(1)]; % Scatterer 3
sc4 = [5; 5; 5+4*rand(1)]; % Scatterer 4

c = 3e8; % speed of light
Q_true = deg2rad(5); % elevation of HV (answer)
w_true = deg2rad(70); % angle of HV compare to SV (answer)

SV = [0; 0; 0]; % Sensing Vehicle position
HVarr = [6,-15,5.5; 8,-15,5.5; 10,-15,5.5; 12,-15,5.5; 14,-15,5.5;
      6,-20,5.1; 8,-20,5.1; 10,-20,5.1; 12,-20,5.1; 14,-20,5.1;
      6,-25,4.7; 8,-25,4.7; 10,-25,4.7; 12,-25,4.7; 14,-25,4.7;
      6,-30,4.3; 8,-30,4.3; 10,-30,4.3; 12,-30,4.3; 14,-30,4.3;
      6,-35,3.9; 8,-35,3.9; 10,-35,3.9; 12,-35,3.9; 14,-35,3.9;
      6,-40,3.5; 8,-40,3.5; 10,-40,3.5; 12,-40,3.5; 14,-40,3.5;
      ];

% Q_state와 w_state 값을 기록할 벡터 초기화
x_estim = zeros(size(HVarr,1),1);
y_estim = zeros(size(HVarr,1),1);
z_estim = zeros(size(HVarr,1),1);

parfor indh = 1:size(HVarr,1)
HV = HVarr(indh,:); % Hidden Vehicle position

% Z angle of AoA
alpha = [atan2(sqrt(sc1(1)^2+sc1(2)^2),sc1(3)); atan2(sqrt(sc2(1)^2+sc2(2)^2),sc2(3)); atan2(sqrt(sc3(1)^2+sc3(2)^2),sc3(3)); atan2(sqrt(sc4(1)^2+sc4(2)^2),sc4(3))];
% AoA
theta = [atan2(sc1(2),sc1(1)); atan2(sc2(2),sc2(1)); atan2(sc3(2),sc3(1)); atan2(sc4(2),sc4(1))];
% Z angle of AoD
psi = [atan2(sqrt((sc1(1)-HV(1))^2+(sc1(2)-HV(2))^2),sc1(3)-HV(3)); atan2(sqrt((sc2(1)-HV(1))^2+(sc2(2)-HV(2))^2),sc2(3)-HV(3)); atan2(sqrt((sc3(1)-HV(1))^2+(sc3(2)-HV(2))^2),sc3(3)-HV(3)); atan2(sqrt((sc4(1)-HV(1))^2+(sc4(2)-HV(2))^2),sc4(3)-HV(3))] - Q_true;
% AoD
phi = [atan2(sc1(2)-HV(2),sc1(1)-HV(1)); atan2(sc2(2)-HV(2),sc2(1)-HV(1)); atan2(sc3(2)-HV(2),sc3(1)-HV(1)); atan2(sc4(2)-HV(2),sc4(1)-HV(1))] - w_true;

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

fprintf('Iteration : %d, Q(t) : %f , w(t) : %f\n',count,rad2deg(Q_state),rad2deg(w_state));

end
elapsed = toc;
avg_iter_time = elapsed / count;
fprintf('계산 소요 시간 : %.6f ms\n ,평균 1회 iteration 시간: %.6f ms\n', elapsed*1000 , avg_iter_time*1000);

x_estim(indh) = sum(v_state.*(sin(alpha).*cos(theta)+sin(psi+Q_state).*cos(phi+w_state))-(d1+c*tdoa).*sin(psi+Q_state).*cos(phi+w_state))/size(v_state,1);
y_estim(indh) = sum(v_state.*(sin(alpha).*sin(theta)+sin(psi+Q_state).*sin(phi+w_state))-(d1+c*tdoa).*sin(psi+Q_state).*sin(phi+w_state))/size(v_state,1);
z_estim(indh) = sum(v_state.*(cos(alpha)+cos(psi+Q_state))-(d1+c*tdoa).*cos(psi+Q_state))/size(v_state,1);

end
figure;

hold on;
grid on;

walls = [
    -5, 5, 0; -5, 5, 18; 5, 5, 18; 5, 5, 0;
    15, 5, 0; 15, 5, 18; 25, 5, 18; 25, 5, 0;
    5, 20, 0; 5, 20, 18; 5, 5, 18; 5, 5, 0;
    15, 20, 0; 15, 20, 18; 15, 5, 18; 15, 5, 0;
    -5, -5, 0; -5, -5, 18; 5, -5, 18; 5, -5, 0;
    15, -5, 0; 15, -5, 18; 15, -40, 18; 15, -40, 0;
    15, -5, 0; 15, -5, 18; 25, -5, 18; 25, -5, 0;
    5, -5, 0; 5, -5, 18; 5, -40, 18; 5, -40, 0
]; % visualize walls

grounds = [
    -5, 5, 0; -5, -5, 0; 18, -5, 0; 18, 5, 0;
    18, -5, 0; 18, 5, 0; 25, 5, 3; 25, -5, 3;
    5, -40, 3.3; 5, -15, 5.3; 15, -15, 5.3; 15, -40, 3.3;
    5, -15, 5.3; 5, -5, 0; 15, -5, 0; 15, -15, 5.3;
    5, -40, 3.3; 5, -40, 0; 5, -5, 0; 5, -15, 5.3;
    15, -40, 3.3; 15, -40, 0; 15, -5, 0; 15, -15, 5.3;
    15, 20, 0; 5, 20, 0; 5, 5, 0; 15, 5, 0
]; % visualize grounds

scatter3(SV(1), SV(2), SV(3), 100, 'r', 'filled', 'DisplayName', 'SV');
scatter3(HVarr(:,1), HVarr(:,2), HVarr(:,3), 100, 'b', 'filled', 'DisplayName', 'HV');

scatter3(x_estim, y_estim, z_estim, 150, 'm', 'filled', 'DisplayName', 'Estimation');

scatter3(sc1(1), sc1(2), sc1(3), 100, 'k', 'filled', 'DisplayName', 'Scatterer 1');
scatter3(sc2(1), sc2(2), sc2(3), 100, 'k', 'filled', 'DisplayName', 'Scatterer 2');
scatter3(sc3(1), sc3(2), sc3(3), 100, 'k', 'filled', 'DisplayName', 'Scatterer 3');
scatter3(sc4(1), sc4(2), sc4(3), 100, 'k', 'filled', 'DisplayName', 'Scatterer 4');

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
legend('SV', 'HV', 'Estimation','Scatterer','FontSize',12);
hold off;


figure;
hold on;
x_err = HVarr(:,1)-x_estim;
y_err = HVarr(:,2)-y_estim;
z_err = HVarr(:,3)-z_estim;

pos_err = sqrt(x_err.^2+y_err.^2+z_err.^2);
scatter3(HVarr(:,1), HVarr(:,2), pos_err, 150, 'm', 'filled', 'DisplayName', 'Estimation');

grounds = [
    5, -40, 3.3; 5, -15, 5.3; 15, -15, 5.3; 15, -40, 3.3;
    5, -15, 5.3; 5, -5, 0; 15, -5, 0; 15, -15, 5.3;
    5, -40, 3.3; 5, -40, 0; 5, -5, 0; 5, -15, 5.3;
    15, -40, 3.3; 15, -40, 0; 15, -5, 0; 15, -15, 5.3
]; % visualize grounds

fill3(grounds(1:4,1), grounds(1:4,2), grounds(1:4,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(grounds(5:8,1), grounds(5:8,2), grounds(5:8,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(grounds(9:12,1), grounds(9:12,2), grounds(9:12,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill3(grounds(13:16,1), grounds(13:16,2), grounds(13:16,3), [0.182,0.126,0.054], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;

figure;
hold on;

pos_err_reshap = reshape(pos_err,5,[]);
pos_err_dist = mean(pos_err_reshap,1)';

plot(flipud(unique(HVarr(:,2))),pos_err_dist);
set(gca, 'XDir', 'reverse');