function [] = build_trans_mat(points_per_box,map_select)
%BUILD_TRANS_MAT Build the transition matrix
%   Save sample points for get_lyap, and the eigenvector

tic;
toc_step=1;
toc_interval=30; %how often the updates show up

disp('Loading recurrent.dat...')
toc

boxes=load('recurrent-16501.dat');%recurrent-rchpm-12900.dat
init_num_boxes=size(boxes);
good_box_index=mode(boxes(:,5));

%delete undesirable boxes based on the index in column 5 of recurrent.dat
bad_boxes=0;
for b=1:init_num_boxes(1)
    if (boxes(b,5)~=good_box_index) %keeping boxes with most common index
        bad_boxes = [bad_boxes,b];
    end
end
bad_boxes(1)=[]; %get rid of the initial dummy 0
boxes(bad_boxes,:)=[];
[num_boxes,num_cols]=size(boxes);

save bad_boxes.dat bad_boxes -ascii

xboxes=boxes(:, [1 2]);
yboxes=boxes(:, [3 4]);
xlength=xboxes(1,2)-xboxes(1,1);
ylength=yboxes(1,2)-yboxes(1,1);

trans_mat=zeros(num_boxes,num_boxes);
sample_points=zeros(num_boxes,2);

disp('Entering loop for building transition matrix... Updates on the minute...')
toc

for i=1:num_boxes
    
%     %set up points to be used for constructing transition matrix
%     xmesh = xboxes(i,1):(xboxes(i,2)-xboxes(i,1))/(sqrt(points_per_box)-1):xboxes(i,2);
%     ymesh = yboxes(i,1):(yboxes(i,2)-yboxes(i,1))/(sqrt(points_per_box)-1):yboxes(i,2);
%     
%     %mesh = zeros(points_per_box,2);
%     counter=1;
%     counter2=1;
    
    for j=1:points_per_box
        %select random points in each box
        rx = xboxes(i,1) + (xboxes(i,2)-xboxes(i,1)).*rand(1,1);
        ry = yboxes(i,1) + (yboxes(i,2)-yboxes(i,1)).*rand(1,1);
        
%         if j<=counter*sqrt(points_per_box)
%             rx = xmesh(counter);
%             ry = ymesh(counter2);
%             counter2=counter2+1;
%         else
%             counter = counter+1;
%             counter2=1;
%             rx = xmesh(counter);
%             ry = ymesh(counter2);
%             counter2=counter2+1;
%         end
        
        
%         %select points in a mesh of the box
%         if mod(j,sqrt(points_per_box))~=0
%             rx = xmesh(mod(j,sqrt(points_per_box)));
%         else
%             rx = xmesh(sqrt(points_per_box));
%         end
%         ry = ymesh(ceil(j/sqrt(points_per_box)));

        %map those random points
        [frx,fry] = the_map(rx,ry,map_select);

        %find out which box they ended up in
        foundx = find(abs(xboxes(:,1)-frx) < xlength & xboxes(:,1)< frx );
        foundy = find(abs(yboxes(:,1)-fry) < ylength & yboxes(:,1)< fry );
        theimage = intersect(foundx, foundy);

        %update the transition matrix
        trans_mat(i, theimage) = trans_mat(i, theimage) + 1;
    end
    
    %print out sample points to be used by get_lyap
    %take center points, or just use last rx,ry
    sample_points(i,:)=[xboxes(i,1)+xlength/2, yboxes(i,1)+ylength/2];
    
    if toc > toc_interval*toc_step
        toc
        toc_step=toc_step+1;
        disp(['Percentage of boxes completed: ', num2str(100*i/num_boxes)]);
    end
end

disp('Done with loop... Calling eigs...')
toc

%make the transition matrix into a probability matrix
trans_mat = trans_mat.*(1/(points_per_box));

disp(['Sum of first row: ', num2str(sum(trans_mat(1,:))) ]);
disp(['Sum of second row: ', num2str(sum(trans_mat(2,:))) ]);

%take transpose of matrix, cast it as sparse,
%then find eigenvector corrsesponding to 1, and normalize it
Final_mat=sparse(trans_mat');

%D might have approximate eigenvalues 1 and -1 
%...but sometimes approx abs(-1) is greater than approx 1 due to roundoff
[V,D]=eigs(Final_mat);
if D(1)>0
    the_eigenvec = V(:,1);
else
    the_eigenvec = V(:,2);
end
the_eigenvec = the_eigenvec./sum(the_eigenvec);

save the_eigenvec.dat the_eigenvec -ascii
save sample_points.dat sample_points -ascii
%save transition_matrix.dat trans_mat -ascii

%the run time
TimeSpent = toc/60;

disp(['Total time to build the transition matrix: ', num2str(TimeSpent), ' minutes.' ]);

end

