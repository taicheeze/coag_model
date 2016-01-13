function coag_plot2(y,species_list,t)
t1=cell2mat(t);
for k=1:length(species_list)
    figure(k)
    plot(t1,y(:,k))
    title(species_list{k})
end