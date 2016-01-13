function [t1, y] = coag_plot(C,t,species_list,dx,dy,cutoff,xrange,x,yrange,y)
C2=C;
area=dx.*dy;
if nargin == 10
    x1=xrange(1);
    x2=xrange(2);
    y1=yrange(1);
    y2=yrange(2);
    [~,indx1]=find(x>=x1,1);
    [~,indx2]=find(x<=x2,1,'last'); 
    [indy1,~]=find(y>=y1,1);
    [indy2,~]=find(y<=y2,1,'last');
    area(:,indx2:end)=[];
    area(:,1:indx1)=[];
    area(indy2:end,:)=[];
    area(1:indy1,:)=[];
    dx(:,indx2:end)=[];
    dx(:,1:indx1)=[];    
    for j=1:cutoff;
        for i=1:length(t);
            C2{i,j}(:,indx2:end)=[];
            C2{i,j}(:,1:indx1)=[];
            C2{i,j}(indy2:end,:)=[];
            C2{i,j}(1:indy1,:)=[];
        end
    end
    for j=cutoff+1:length(species_list);
        for i=1:length(t);
            C2{i,j}(:,indx2:end)=[];
            C2{i,j}(:,1:indx1)=[];
        end
    end
end
totarea=sum(sum(area));
totlength=dx(1)*size(area,2);
t1=cell2mat(t);
y=zeros(size(C2));
for j=1:cutoff;
    for i=1:length(t);
        amount=sum(sum(C2{i,j}.*area));
        y(i,j)=amount/totarea;
    end
end
for j=cutoff+1:length(species_list);
    for i=1:length(t);
        amount=sum(C2{i,j}.*dx(end,:));
        y(i,j)=amount/totlength;
    end
end
for k=1:length(species_list)
    figure(k)
    plot(t1,y(:,k))
    title(species_list{k})
end
        
% Z=[t1, 1e6*(y(:,2)+1.2*y(:,5)), 1e6*y(:,3), 1e6*y(:,1)];
 %%
% csvwrite('PDE_IIa_FVa_control.csv',Z);