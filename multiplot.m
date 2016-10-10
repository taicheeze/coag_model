function multiplot(col,t,s,A,B,C,D);
plot(t,A(:,col),t,B(:,col),t,C(:,col),t,D(:,col));
title(s(col));
legend('plot1','plot2','plot3','plot4')