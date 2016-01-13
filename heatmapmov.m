function heatmapmov(C,X,Y,t,filename)
close all hidden
C1=cell2mat(C);
C1=reshape(C1,size(X,1),size(X,2),length(t));
cmax=max(max(max(C1)));
for i=2:length(t)
    contourf(X,Y,C{i})
    caxis([0 cmax]);
    colorbar
    time=num2str(t{i});
    title(['time= ' time]);
    set(gca,'ydir','normal');
    %drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i==2;
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    close all hidden
end