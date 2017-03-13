%Calculate measure entropy

%... Make changes in this block...
map_select=6; %see the_map.m - 1,2,3,4,5,6
points_per_box=100; %used in building the probability transition matrix
                    % ... use a square amount for a mesh instead of random
lyap_method=4; %see get_lyap.m - 1,2,3,4 and 5 is a little different
delta=1.0e-11; %used for lyap_method=2
max_iters=1;%10000; %lyap iterations, set to 1 for local lyap if vector is ok
do_build_mat=true; %call build_trans_mat.m or not
%... everything else should take care of itself...
%... unless you change the file names that are loaded...

%transition matrix stuff happens here if desired
if(do_build_mat)
    build_trans_mat(points_per_box,map_select);
end

%start timer
tic;

%that just gave us the eigenvector and the sample points for get_lyap
the_eigenvec=load('the_eigenvec.dat'); %the_eigenvec_large_henon.dat
sample_points=load('sample_points.dat');%sample_points_large_henon.dat
num_boxes=length(the_eigenvec);
lyapexp=zeros(num_boxes,1);

toc_step=1;
toc_interval=60; %how often the updates show up

disp('Entering loop for finding exponents... Updates on the minute...');
disp('If it appears to be stopped, then theres a blow up in the function or roundoff.');

for i=1:num_boxes
    lyapexp(i) =... 
        get_lyap(sample_points(i,1),sample_points(i,2),delta,map_select,lyap_method,max_iters);
    
    if toc > toc_interval*toc_step
        toc
        toc_step=toc_step+1;
        disp(['Percentage of boxes completed: ', num2str(100*i/num_boxes)]);
    end
end

disp('Done with loop...');

%the entropy
entropy = sum(the_eigenvec.*lyapexp);

%the run time
TimeSpent = toc;

disp(['Total time (not including building the transition matrix): ',...
    num2str(TimeSpent), ' seconds.' ]);
disp(['Minimum exponent: ', num2str(min(lyapexp)) ]);
disp(['Maximum exponent: ', num2str(max(lyapexp)) ]);
disp(['Entropy: ', num2str(entropy) ]);
