%Plot invariant measure

the_eigenvec=load('the_eigenvec.dat'); 
boxes=load('recurrent-16501.dat');
bad_boxes=load('bad_boxes.dat');

boxes(bad_boxes,:)=[];
X=boxes(:,[1,2]);
Y=boxes(:,[3,4]);
Z=the_eigenvec;
Z=abs(Z);

p=patch([X(:,1),X(:,2),X(:,2),X(:,1)]',[Y(:,1),Y(:,1),Y(:,2),Y(:,2)]',[Z,Z,Z,Z]','red')
set(p,'EdgeColor','b')
%set(p,'EdgeAlpha',0.1)
%set(p,'FaceColor','none')
