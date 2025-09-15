%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Title  : Iterative NLoS 투시 알고리즘에서 코드 단순화를 위한 삼각함수 연산 과정 util 함수
% Type   : Function
% Input  : opt, val1, val2
% Output : opt의 값에 따른 삼각함수 연산 네가지 옵션 중 택 1 실행한 결과값

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = triop(opt,val1,val2)

if (opt == 1)
val = sin(val1).*cos(val2);
elseif (opt == 2)
val = sin(val1).*sin(val2);
elseif (opt == 3)
val = cos(val1).*cos(val2);
elseif (opt == 4)
val = cos(val1).*sin(val2);
else
fprintf('Incorrect value');
pause;
end
end