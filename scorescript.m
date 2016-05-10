normvals=zeros(2,5);
stdvals=zeros(1,5);
for i=1:5
    normvals(1,i)=norm(scores(:,i)-scores(:,6),1);
    normvals(2,i)=norm(scores(:,i)-scores(:,6),2);
    stdvals(1,i)=std(abs(scores(:,i)-scores(:,6)));
end
differences=zeros(30,5);
for j=1:5
    differences(:,j)=abs(scores(:,j)-scores(:,6));
end
teamstd=zeros(30,1);
teamavg=zeros(30,1);
teamrange=zeros(30,1);
for i=1:30
    teamstd(i)=std(differences(i,:));
    teamavg(i)=mean(differences(i,:));
    teamrange(i)=range(scores(i,1:5));
end
    
