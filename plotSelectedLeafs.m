function toplot = plotSelectedLeafs(selectedLeafs, runPath, runName, MRoptions)
% currently plots the bar graph of MMPConReKS

% specify the colors of :(< threshold, > threshold, NaN) in bar graph
colors = {[0 154 205]/255;[139 28 98]/255;[238 238 0]/255};

toplot = [];
for i=1:length(selectedLeafs.ID)
    toplot(i,:) = [sum(selectedLeafs.flag{i}==1) sum(selectedLeafs.flag{i}==0) sum(isnan(selectedLeafs.flag{i}))]; 
end

hasNaN = logical(sum(toplot(:,3)));

figure
if hasNaN
    hBar = bar(toplot,'stacked')
    set(hBar,{'FaceColor'},colors);
    legend(['ind(x";T|[cs\x] p < ' num2str(MRoptions.threshold1) ' & ind(x;T|x" p > ' num2str(MRoptions.threshold2)],['ind(x";T|[cs\x] p >= ' num2str(MRoptions.threshold1)  ' & ind(x;T|x" p <= ' num2str(MRoptions.threshold2)], 'NaN')
else
    toplot(:,3)=[];
    hBar = bar(toplot,'stacked')
    set(hBar,{'FaceColor'},colors(1:2));
    legend(['ind(x";T|[cs\x] p < ' num2str(MRoptions.threshold1) ' & ind(x;T|x" p > ' num2str(MRoptions.threshold2)],['ind(x";T|[cs\x] p >= ' num2str(MRoptions.threshold1)  ' & ind(x;T|x" p <= ' num2str(MRoptions.threshold2)])
end

title(runName)
xlabel('selected variables')
ylabel('level counting from leaf')
set(gca,'XTickLabel',selectedLeafs.varNames(selectedLeafs.ID)')
saveas(gcf,[runPath runName],'fig')
