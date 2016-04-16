function multiplot(col,t,s,A,B,C,D,E);
plot(t,A(:,col),t,B(:,col),t,C(:,col),t,D(:,col),t,E(:,col));
title(s(col));
legend('plot1','plot2','plot3','plot4','plot5')