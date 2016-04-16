function [C2, X, Y]=matresize(C,col,x,y,x1,x2,y1,y2)
if nargin==6||8
    C2=C(:,col);
    X=x;
    Y=y;
    [~,indx1]=find(x>=x1,1);
    [~,indx2]=find(x<=x2,1,'last'); 
    X(:,indx2+1:end)=[];
    X(:,1:indx1-1)=[];
    Y(:,indx2+1:end)=[];
    Y(:,1:indx1-1)=[];
        for i=1:length(C)
            C2{i}(:,indx2+1:end)=[];
            C2{i}(:,1:indx1-1)=[];
        end
    if nargin==8;
        [indy1,~]=find(y>=y1,1);
        [indy2,~]=find(y<=y2,1,'last');
        X(indy2+1:end,:)=[];
        X(1:indy1-1,:)=[];
        Y(indy2+1:end,:)=[];
        Y(1:indy1-1,:)=[];
        for i=1:length(C)
            C2{i}(indy2+1:end,:)=[];
            C2{i}(1:indy1-1,:)=[];
        end
    end
else
    error('wrong number of inputs')
end