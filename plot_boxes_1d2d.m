function z=plot_boxes_1d2d(file,ax,col,fig,stretch)
% file is the filename in single quotes 'filename.dat' to be loaded.
% ax is a vector for the plotting axis in the form [xmin xmax ymin ymax]
% col is the color of the boxes such as 'b' , 'r' , 'm' , 'k' etc.
% fig is the number of the figure in which to draw the boxes
% stretch is the (optional) stretch factor for the height of the boxes when plotting in 1d 

if nargin==4
    stretch=1;
end

figure(fig);
axis(ax);
hold on;

boxes=load(file,'-ascii');
[sb1,sb2]=size(boxes);

hp=patch([boxes(:,1) boxes(:,2) boxes(:,2) boxes(:,1)]',stretch*[boxes(:,3) boxes(:,3) boxes(:,4) boxes(:,4)]',col);
set(hp,'EdgeColor','none');
