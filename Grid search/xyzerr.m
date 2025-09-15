function err = xyzerr(answer,estimM)
errM = estimM-answer;

err = figure;
subplot(1,3,1)
histogram(errM(:,1),'BinWidth', 0.5)
title('Error histogram of x');
subplot(1,3,2)
histogram(errM(:,2),'BinWidth', 0.5)
title('Error histogram of y');
subplot(1,3,3)
histogram(errM(:,3),'BinWidth', 0.5)
title('Error histogram of z');

end
