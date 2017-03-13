function [u,v] = the_map(x,y,selection)
%THE_MAP Iterate the map
%   1 - reinjected horseshoe map
%   2 - reinjected cuspidal horseshoe Poincare map
%   3 - Henon map
%   4 - g2 1D map

if(selection==1) %reinjected horseshoe map
    %.59,.2,1.3 looks like a long periodic orbit
    cusp = 0.65;%1.5;%.95;%0.65;
    u2 = 0;
    u3 = 1/2;
    u4 = 1;
    u5 = 3/2;
    u6 = 0.2;%;%.2;
    u7 = 0.6;%1.5;%.6;%.6;
    v0 = 0;
    v1 = 1/3;
    v2 = 2/3;
    v3 = 1;
    v4 = .4;
    v5 = .6;
    midv = (v3+v0)/2;
    radv = midv-(v1-v0)/v3.*y;
    lefv = midv + radv;
    rigv = midv - radv;
    u1 = -midv;

    if u1<=x && x<u2 %LHS
        u = (u6-u7)/(u1-u2).*(x-u1)+u6;
        v = (v4-v5)/(v0-v3).*(y-v0)+v4;

    elseif u2<=x && x<=u3 %LRHS
        u = (cusp-u2)/(u2-u3).*(x-u2)+cusp;
        v = (v3-v2)/(v0-v3).*(y-v0)+v3;

    elseif u3<x && x<u4 %MRHS
        %Original Calculation
        %v = (lefv-rigv)/(u3-u4).*(x-u3)+lefv; %used for u calculation
        %u = -( (radv).^2 - (v-midv).^2 ).^(1/2);
        
        v =((4*y)/3 - 2)*(x - 1/2) - y/3 + 1;
        % expanded (4*x*y)/3 - y - 2*x + 2
        
        u = -((y/3 - 1/2)^2 - (((4*y)/3 - 2)*(x - 1/2) - y/3 + 1/2)^2)^(1/2);
        %expanded -(- (16*x^2*y^2)/9 + (16*x^2*y)/3 - 4*x^2 + (8*x*y^2)/3 - 8*x*y + 6*x - (8*y^2)/9 + (8*y)/3 - 2)^(1/2)

    elseif u4<=x && x<=u5  %RRHS
        u = (u2-u5)/(u4-u5).*(x-u4)+u2;
        v = (v0-v1)/(v0-v3).*(y-v0);

    else
        disp('OUT OF BOUNDS')
    end

elseif(selection==2) %reinjected cuspidal horseshoe poincare map
    pm_c=2.75;
    pm_M=3.47098748700613;
    pm_epsilon=0.025;
    pm_alpha=2.1;
    pm_a1=6;
    pm_w1=0.3;
    pm_k1=-3.0;
    pm_d1=-2.0;
    pm_b=0.4;
    pm_a2=1560.0;
    pm_c2=0.4;
    pm_w2=-0.325;
    pm_k2=6.0;
    pm_d2=-2.0;
    pm_E=1.7; %// used to be 1.0
    pm_plus=0.25;

    %DOMAIN : [-0.13125,2.66875] X [-0.48,0.31]

    pm_eps_alpha=pm_epsilon^pm_alpha;%exp(pm_alpha*ln(pm_epsilon));
    pm_eps_alpha1=pm_epsilon^(pm_alpha-1);%exp((pm_alpha-one)*ln(pm_epsilon));

    if(x<=0)
        u=pm_E-pm_b*sqrt(abs(x))/sqrt(pm_epsilon)-pm_c2*x+pm_a2*abs(x)^pm_alpha*y;
        v=pm_w2+(pm_k2*sqrt(pm_epsilon)+pm_plus)*sqrt(abs(x))+pm_d2...
            *pm_eps_alpha*abs(x)^pm_alpha*y;
    elseif(x>0)
        u=1-pm_M*sqrt(x)+pm_c*x-pm_a1*pm_eps_alpha1*x^pm_alpha*y;
        v=pm_w1+pm_k1*sqrt(pm_epsilon)*sqrt(x)+pm_d1*pm_eps_alpha*x^pm_alpha*y;
    else
        disp('BAD RCH_PM EVALUATION');
    end
    
elseif(selection==3) %Henon map
    u=1-1.4.*x.^2+y;
    v=0.3.*x;
    
elseif(selection==4) %g2 map
    if(x<=1/3)
        u=2*x/(1-x);
        v=0;
    elseif(x<=1)
        u=(1-x)/(2*x);
        v=0;
    else
        disp('OUT OF BOUNDS');
    end
    
elseif(selection==5)
    cusp = 0.65;%1.5;%0.65;%
    u2 = 0;
    u3 = 1/2;
    u4 = 1;
    u5 = 3/2;
    u6 = 0.2;%0;%.2;%
    u7 = 0.6;%1.5;%.6;%
    v0 = 0;
    v1 = 1/3;
    v2 = 2/3;
    v3 = 1;
    v4 = .4;
    v5 = .6;
    midv = (v3+v0)/2;
    u1 = -midv; %has to be the max radius of semicircle
    
    if u1<=x && x<u2 %LHS
        u = (u6-u7)/(u1-u2).*(x-u1)+u6;
        v = (v4-v5)/(v0-v3).*(y-v0)+v4;
        
    elseif u2<=x && x<=u3 %LRHS
        u = (cusp-u2)/(u2-u3).*(x-u2)+cusp;
        v = (v3-v2)/(v0-v3).*(y-v0)+v3;
        
    elseif u3<x && x<u4 %MRHS
        %Original Calculation
        %r = -1/3*y+1/2;
        %theta = (pi/2-3/2*pi)/(u3-u4)*(u-u3)+pi/2;
        %u = r*cos(theta);
        %v = r*sin(theta)+midv;
        
        u = -cos(pi/2 + 2*pi*(x - 1/2))*(y/3 - 1/2);
        v = 1/2 - sin(pi/2 + 2*pi*(x - 1/2))*(y/3 - 1/2);
        
    elseif u4<=x && x<=u5  %RRHS
        u = (u2-u5)/(u4-u5).*(x-u4)+u2;
        v = (v0-v1)/(v0-v3).*(y-v0);
        
    else
        disp('OUT OF BOUNDS')
    end
    
elseif(selection==6) %New rchpm map
    c=2.75;
    M=3.47098748700613;
    eps=0.025;
    alph=2.1;
    
    if(x>0)
        u= 1-M*sqrt(x)+c*x-6*eps^(alph-1)*(x^alph)*y;
        v= 0.3-3*sqrt(eps*x)-2*(eps^alph)*(x^alph)*y;
    elseif(x<0)
        %original
        %u= 0.65+0.4*sqrt(abs(x)/eps)+0.4*x+240*(abs(x)^alph)*y;
        %v= -0.225+(6*sqrt(eps)+0.25)*sqrt(abs(x))-2*(eps^alph)*(abs(x)^alph)*y;
        
        u= 0.65+0.4*sqrt(-x/eps)+0.4*x+240*((-x)^alph)*y;
        v= -0.225+(6*sqrt(eps)+0.25)*sqrt(-x)-2*(eps^alph)*((-x)^alph)*y;
    else
        disp('BAD RCH_PM EVALUATION');
        
    end
    
end %end selection if statements

end
