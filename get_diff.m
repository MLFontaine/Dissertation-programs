function [ J ] = get_diff(x,y,map_select)
%GET_DIFF For use in get_lyap
%   Need the parameters to match up with the_map.m

if(map_select==1)
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
        
        J = [(u6-u7)/(u1-u2), 0; 0, (v4-v5)/(v0-v3)];
        
    elseif u2<=x && x<=u3 %LRHS
        
        J = [(cusp-u2)/(u2-u3), 0; 0, (v3-v2)/(v0-v3)];
        
    elseif u3<x && x<u4 %MRHS
        
        J = [(((4*y)/3 - 2)*(((4*y)/3 - 2)*(x - 1/2) - y/3 + 1/2))/((y/3 - 1/2)^2 - (((4*y)/3 - 2)*(x - 1/2) - y/3 + 1/2)^2)^(1/2), (2*((4*x)/3 - 1)*(((4*y)/3 - 2)*(x - 1/2) - y/3 + 1/2) - (2*y)/9 + 1/3)/(2*((y/3 - 1/2)^2 - (((4*y)/3 - 2)*(x - 1/2) - y/3 + 1/2)^2)^(1/2));
            (4*y)/3 - 2, (4*x)/3 - 1];
        
    elseif u4<=x && x<=u5  %RRHS
        
        J = [(u2-u5)/(u4-u5), 0; 0, (v0-v1)/(v0-v3)];
        
    else
        disp('OUT OF BOUNDS')
    end
    
elseif(map_select==3)
    J = [-2.8.*x, 1; 0.3, 0];
    
elseif(map_select==5)
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
        
        J = [(u6-u7)/(u1-u2), 0; 0, (v4-v5)/(v0-v3)];
        
    elseif u2<=x && x<=u3 %LRHS
        
        J = [(cusp-u2)/(u2-u3), 0; 0, (v3-v2)/(v0-v3)];
        
    elseif u3<x && x<u4 %MRHS
        
        J = [2*pi*sin(pi/2 + 2*pi*(x - 1/2))*(y/3 - 1/2), -cos(pi/2 + 2*pi*(x - 1/2))/3;
             -2*pi*cos(pi/2 + 2*pi*(x - 1/2))*(y/3 - 1/2), -sin(pi/2 + 2*pi*(x - 1/2))/3];
        
        %J = [-sin(2*pi*(x - 1/2) - (3*pi)/2)*(1/2 - y/3)*2*pi, cos(2*pi*(x - 1/2) - (3*pi)/2)*(-1/3);
        %    -cos(2*pi*(x - 1/2) - (3*pi)/2)*(1/2 - y/3)*2*pi, -sin(2*pi*(x - 1/2) - (3*pi)/2)*(-1/3)];
        
    elseif u4<=x && x<=u5  %RRHS
        
        J = [(u2-u5)/(u4-u5), 0; 0, (v0-v1)/(v0-v3)];
        
    else
        disp('OUT OF BOUNDS')
    end
    
    
elseif(map_select==6)
    c=2.75;
    M=3.47098748700613;
    eps=0.025;
    alph=2.1;
    
    if(x>=0)
        J = [c - M/(2*x^(1/2)) - 6*alph*eps^(alph - 1)*x^(alph - 1)*y, -6*eps^(alph - 1)*x^alph;
             - (3*eps)/(2*(eps*x)^(1/2)) - 2*alph*eps^alph*x^(alph - 1)*y, -2*eps^alph*x^alph];
        
        %J = [ -M*0.5*(1/sqrt(x)) + c - 6*eps^(alph-1)*alph*(x^(alph-1))*y, -6*eps^(alph-1)*(x^alph);
        %    -3*0.5*eps*(1/sqrt(eps*x)) - 2*(eps^alph)*alph*(x^(alph-1))*y, -2*(eps^alph)*(x^alph)];
        
    elseif(x<0)
        J = [2/5 - 240*alph*(-x)^(alph - 1)*y - 1/(5*eps*(-x/eps)^(1/2)), 240*(-x)^alph;
             2*alph*eps^alph*(-x)^(alph - 1)*y - (6*eps^(1/2) + 1/4)/(2*(-x)^(1/2)), -2*eps^alph*(-x)^alph];
          
        %J= [ 0.4*(-0.5/eps)*(1/sqrt(abs(x)/eps)) + 0.4 + 240*(-1*alph)*(abs(x)^(alph-1))*y, 240*(abs(x)^alph);
        %    (6*sqrt(eps)+0.25)*(-0.5)*(1/sqrt(abs(x))) - 2*(eps^alph)*(-1*alph)*(abs(x)^(alph-1))*y, -2*(eps^alph)*(abs(x)^alph)];
        
    end
    
    
end

end

