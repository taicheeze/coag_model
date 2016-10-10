function VX_plot(t,y)
VXtot=y(:,18)+y(:,15)+y(:,16)+y(:,17);
plot(t,VXtot);
title('Prothrombinase')
xlabel('time (sec)')
ylabel('Concentration (mol/m^2)')
figure
plot(t,y(:,12))
title('TF_VIIa')
xlabel('time (sec)')
ylabel('Concentration (mol/m^2)')
figure
plot(t,y(:,10))
title('Xa')
xlabel('time (sec)')
ylabel('Concentration (mol/m^3)')
figure
plot(t,y(:,8))
title('Thrombin')
xlabel('time (sec)')
ylabel('Concentration (mol/m^3)')
figure
plot(t,y(:,3))
title('Fibrinogen')
xlabel('time (sec)')
ylabel('Concentration (mol/m^3)')
figure
plot(t,y(:,4))
title('Fibrin')
xlabel('time (sec)')
ylabel('Concentration (mol/m^3)')

figure
plot(t,y(:,8)/max(y(:,8)),t,VXtot/max(VXtot),t,y(:,12)/max(y(:,12)),t,y(:,10)/max(y(:,10)),t,y(:,9)/max(y(:,9)))
legend('thrombin','prothrombinase','TF_VIIa','Xa','mIIa')