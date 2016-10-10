lengths=[10 25 50 100];
injs=[6 5.4 6.4 7.0];
%k=12;
for i=1:3%:length(injs)
    for j=1:length(lengths)
        name=['size' num2str(lengths(j)) '_inj_' num2str(i)];
        %eval(['[t, ' name ', s, x, y, dx, dy]=custom_solver_data_analysis(300,' num2str(k) ');']);
        name2=[name '_1'];
        eval(['[t1, ' name2 ']=coag_plot(' name ',t,s,dx,dy,9,[5e-3 6e-3],x,[0 10e-6],y);']);
        %disp(['parsed folder ' num2str(k)])
        %k=k+1;
    end
end