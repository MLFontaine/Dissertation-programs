function [lyap] = get_lyap(x,y,delta,map_select,method_select,max_iters)
%GET_LYAP Get the Lyapunov exponent for the given function

if(method_select==1) %Call TISEAN
    tiseanPath = '/Users/MF/bin/';
    
    s(1,1) = x; s(1,2) = y;
    for i = 2:10001 %plot out 10,000 points... more points seems to be better
        [s(i,1),s(i,2)] = the_map(s(i-1,1),s(i-1,2),map_select);
    end
    s(1,:) = [];
    save the_map.dat s -ascii -tabs
    %plot(s(:,1),s(:,2),'.')
    
    %need to have TISEAN installed for this system call
    system([tiseanPath,'lyap_r -V0 -s10 -o lyap.dat the_map.dat']);
    l = load('lyap.dat');
    lyapfit = polyfit(l(2:10,1), l(2:10,2),1);
    lyap = lyapfit(1);
    
elseif(method_select==2) %One separation
    nbrhd_type=2; %pick the type of delta nbrhd
    %%%max_iters=100000;%1/(delta); %10000
    
    if(nbrhd_type==1) %take points inside delta circle
        xx=100;
        yy=100;
        while(sqrt(abs(x-xx)^2+abs(y-yy)^2)>delta)
            xx=x - delta/2 + delta*rand;
            yy=y - delta/2 + delta*rand;
        end
    elseif(nbrhd_type==2) %one of 4 points exactly delta away
        randx=(rand-1/2);
        randy=(rand-1/2);
        xx = x+randx/abs(randx)*delta/sqrt(2);
        yy = y+randy/abs(randy)*delta/sqrt(2);
    elseif(nbrhd_type==3) %take points inside square with sides of delta
        xx = x - delta/2 + delta*rand;
        yy = y - delta/2 + delta*rand;
    end
    
    for i=1:max_iters
        [x,y] = the_map(x,y,map_select);
        [xx,yy] = the_map(xx,yy,map_select);
        if( sqrt(abs(x-xx)^2+abs(y-yy)^2)>=1 )
            break;
        end
    end
    
    if(i<max_iters)
        lyap = -log(delta)/i;
    else
        lyap=NaN;
        disp('Failure to separate.')
    end
    
elseif(method_select==3) %Exact calculation
    %%%max_iters=10000;
    store = zeros(1,max_iters);
    mysave = zeros(1,max_iters);
    my_det = zeros(1,max_iters);
    %pick an arbitrary vector to begin with,
    %might want to change with map_select
    pre_vec = [1;2];
    pre_vec = pre_vec/norm(pre_vec);
    [x,y] = the_map(x,y,map_select);
    for i = 1:max_iters
        J = get_diff(x,y,map_select);
        post_vec = J*pre_vec;
        store(i) = norm(post_vec);
        pre_vec = post_vec/norm(post_vec);
        [x,y] = the_map(x,y,map_select);
        
        mysave(i) = mean(log(store(1:i)));
        my_det(i) = log(abs(det(J)));
    end
    lyap = mean(log(store));
    
    disp(['lyap2=  ', num2str(mean(my_det)-lyap)]);
    plot_num=max_iters-1;
    plot((max_iters-plot_num:max_iters),mysave((max_iters-plot_num):max_iters),'-');
    
elseif(method_select==4) %Max out the direction vector in one iteration
    %find the max expansion direction
    angs = 0:.001:2*pi;
    store = zeros(1,length(angs));
    for k = length(angs)
        pre_vec = [cos(angs(k));sin(angs(k))];
        pre_vec = pre_vec/norm(pre_vec);
        [x,y] = the_map(x,y,map_select);
        J = get_diff(x,y,map_select);
        post_vec = J*pre_vec;
        store(k) = norm(post_vec);
    end
    lyap = max(log(store));
    
elseif(method_select==5) %Compute all lyaps with QR
    
    J = get_diff(x,y,map_select);
    [Q,R]=qr(J);
    Rsave = R;
    for i = 1:max_iters
        [x,y] = the_map(x,y,map_select);
        J = get_diff(x,y,map_select);
        JJ = J*Q;
        [Q,R] = qr(JJ);
        Rsave = R*Rsave;
    end
    lyap1 = log(abs(Rsave(1,1)))/max_iters;
    lyap2 = log(abs(Rsave(2,2)))/max_iters;
    lyap= lyap1+lyap2;
    
    disp(['lyap1=  ', num2str(lyap1)]);
    disp(['lyap2=  ', num2str(lyap2)]);
    
else
    disp('Bad lyap method selection.')
end

end




