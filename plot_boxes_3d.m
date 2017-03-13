function z=plot_boxes_3d(file,ax,col,fig,edge_alpha)
% file is the filename in single quotes 'filename.dat' to be loaded.
% ax is a vector for the plotting axis in the form [xmin xmax ymin ymax]
% col is the color of the boxes such as 'b' , 'r' , 'm' , 'k' etc.
% fig is the number of the figure in which to draw the boxes
% edge_alpha determines the opacity of the edge lines -- default is 0.1

if (nargin==4)
	edge_alpha=0.1;
end

figure(fig);
axis(ax);
hold on;

boxes=load(file,'-ascii');
[sb1,sb2]=size(boxes);

p=patch([boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1) boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1)]',[boxes(:,3) boxes(:,3) boxes(:,4) boxes(:,4) boxes(:,3) boxes(:,3) boxes(:,4) boxes(:,4)]',[boxes(:,5) boxes(:,5) boxes(:,5) boxes(:,5) boxes(:,5) boxes(:,5) boxes(:,5) boxes(:,5)]',col)
set(p,'Edgecolor',col)
set(p,'EdgeAlpha',edge_alpha)
set(p,'FaceColor','none')

p=patch([boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1) boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1)]',[boxes(:,3) boxes(:,3) boxes(:,4) boxes(:,4) boxes(:,3) boxes(:,3) boxes(:,4) boxes(:,4)]',[boxes(:,6) boxes(:,6) boxes(:,6) boxes(:,6) boxes(:,6) boxes(:,6) boxes(:,6) boxes(:,6)]',col)
set(p,'EdgeColor',col)
set(p,'EdgeAlpha',edge_alpha)
set(p,'FaceColor','none')

p=patch([boxes(:,1) boxes(:,1) boxes(:,1) boxes(:,1) boxes(:,1) boxes(:,1) boxes(:,1) boxes(:,1)]',[boxes(:,3) boxes(:,4) boxes(:,4) boxes(:,3) boxes(:,3) boxes(:,4) boxes(:,4) boxes(:,3)]',[boxes(:,5) boxes(:,5) boxes(:,6) boxes(:,6) boxes(:,5) boxes(:,5) boxes(:,6) boxes(:,6)]',col)
set(p,'EdgeColor',col)
set(p,'EdgeAlpha',edge_alpha)
set(p,'FaceColor','none')

p=patch([boxes(:,2) boxes(:,2) boxes(:,2) boxes(:,2) boxes(:,2) boxes(:,2) boxes(:,2) boxes(:,2)]',[boxes(:,3) boxes(:,4) boxes(:,4) boxes(:,3) boxes(:,3) boxes(:,4) boxes(:,4) boxes(:,3)]',[boxes(:,5) boxes(:,5) boxes(:,6) boxes(:,6) boxes(:,5) boxes(:,5) boxes(:,6) boxes(:,6)]',col)
set(p,'EdgeColor',col)
set(p,'EdgeAlpha',edge_alpha)
set(p,'FaceColor','none')

p=patch([boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1) boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1)]',[boxes(:,3) boxes(:,3) boxes(:,3) boxes(:,3) boxes(:,3) boxes(:,3) boxes(:,3) boxes(:,3)]',[boxes(:,5) boxes(:,5) boxes(:,6) boxes(:,6) boxes(:,5) boxes(:,5) boxes(:,6) boxes(:,6)]',col)
set(p,'EdgeColor',col)
set(p,'EdgeAlpha',edge_alpha)
set(p,'FaceColor','none')

p=patch([boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1) boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1)]',[boxes(:,4) boxes(:,4) boxes(:,4) boxes(:,4) boxes(:,4) boxes(:,4) boxes(:,4) boxes(:,4)]',[boxes(:,5) boxes(:,5) boxes(:,6) boxes(:,6) boxes(:,5) boxes(:,5) boxes(:,6) boxes(:,6)]',col)
set(p,'EdgeColor',col)
set(p,'EdgeAlpha',edge_alpha)
set(p,'FaceColor','none')

