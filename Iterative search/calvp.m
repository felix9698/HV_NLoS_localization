%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Title  : Iterative NLoS 투시 알고리즘에서 Q(t)와 w(t)를 기반으로 vp(t) 업데이트 하는 함수
% Type   : Function
% Input  : t번째 iter에서의 Q,W 값, 주어진 P*1 행렬들 (alpha,theta,psi,phi,tdoa), d1 값
% Output : t번째 iter에서의 P*1 v행렬 (v1,...,vp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vt = calvp(Qt,wt,alpha,theta,psi,phi,tdoa,d1,c)

%%% matrix A %%%
A = matA_calvp(alpha, theta, psi, phi, wt, Qt);
%%% matrix B %%%
B = matB_calvp(c, tdoa, psi, Qt, wt, phi, d1);

vt = inv(A.'*A)*A.'*B; % Psuedo-Inverse

end