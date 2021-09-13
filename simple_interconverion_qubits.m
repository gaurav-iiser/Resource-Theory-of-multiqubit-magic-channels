%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Interconversion among qubits%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define (r1,r2,r3), (s1,s2,s3) such that these vectors belong to Bloch
%%% sphere and r1 <= r2, r3, and s1 <= s2,s3.

clear all
format short
t=0;
c=0;
b=0;
p = 0.83
for k = 1:100
    p1= rand()*(1/3);
    p2= rand()*(2/3);
    p3 = 1 - p1 - p2;

    v = sort([sqrt(p1),sqrt(p2),sqrt(p3)]);

    r1 = 0;%v(1);%1/sqrt(6);
    r2 = 1/sqrt(2);%v(3);%1/sqrt(6);
    r3 = 1/sqrt(2);%v(2);%sqrt(2)/sqrt(3);
   

    p1= rand()*(1/3);
    p2= rand()*(2/3);
    p3 = 1 - p1 - p2;

    v = sort([sqrt(p1),sqrt(p2),sqrt(p3)]);
    %p=rand();
    
    s1 = 0;%p*v(1);%0.9/sqrt(3);
    s2 = 0;%p*v(3);%1/2;%sqrt(5)/sqrt(6); 
    s3 = 1;%p*v(2);%(p*sqrt(2))/sqrt(3);%0.9/sqrt(3);%sqrt(5)/sqrt(6);


    phi(:,1) = [1;0;0];
    phi(:,2) = [0;1;0];
    phi(:,3) = [0;0;1];

    rho(:,1) = [r1; r2; r3]; %%% We're considering rho(:,1) to be in the X-subset of the %%%positive octant %%%
    rho(:,2) = [r3;r1;r2];
    rho(:,3) = [r2; r3; r1];
    rho(:,4) = [-r1;r3;r2];
    rho(:,5) = [-r2;r1;r3];
    rho(:,6) = [r2;-r1;r3];
    rho(:,7) = phi(:,3);
    rho(:,8) = phi(:,2);
    rho(:,9) = [-r3;r2;r1];
    rho(:,10)= [r3;r2;-r1];
    
    %rho(:,5) = phi(:,2);
    %rho(:,6) = phi(:,3);

    sigma(:,1) = [s1;s2;s3]; %%% We're considering sigma(:,1) to be in the X-subset of the %%%positive octant %%%
    sigma(:,2) = [s2;s3;s1];
    sigma(:,3) = [s3;s1;s2];
    sigma(:,4) = [-s1;s3;s2];
    sigma(:,5) = phi(:,2);
    sigma(:,6) = phi(:,3);

    %{ 
    %% old %%%
    [mid_vec1, dec1] =  perp_vec(rho(:,1),rho(:,2),rho(:,6));
    [mid_vec2, dec2] =  perp_vec(rho(:,1),rho(:,2),rho(:,3));
    [mid_vec3, dec3] =  perp_vec(rho(:,1),rho(:,4),rho(:,6));
    [mid_vec4, dec4] =  perp_vec(rho(:,1),rho(:,3),rho(:,5));
    [mid_vec5, dec5] =  perp_vec(rho(:,1),rho(:,4),rho(:,5));
    %}
    
    [mid_vec1, dec1] =  perp_vec(rho(:,1),rho(:,6),rho(:,7));
    [mid_vec2, dec2] =  perp_vec(rho(:,1),rho(:,7),rho(:,5));
    [mid_vec3, dec3] =  perp_vec(rho(:,1),rho(:,5),rho(:,4));
    [mid_vec4, dec4] =  perp_vec(rho(:,1),rho(:,4),rho(:,3));
    [mid_vec5, dec5] =  perp_vec(rho(:,1),rho(:,3),rho(:,2));
    [mid_vec6, dec6] =  perp_vec(rho(:,1),rho(:,2),rho(:,6));
    [mid_vec7, dec7] =  perp_vec(rho(:,3),rho(:,4),rho(:,8));
    
    
    [mid_vec8, dec8] =  perp_vec(rho(:,1),rho(:,3),rho(:,2));
    [mid_vec9, dec9] =  perp_vec(rho(:,1),rho(:,2),rho(:,7));
    [mid_vec10, dec10] =  perp_vec(rho(:,1),rho(:,7),rho(:,4));
    [mid_vec11, dec11] =  perp_vec(rho(:,1),rho(:,4),rho(:,8));
    [mid_vec12, dec12] =  perp_vec(rho(:,1),rho(:,8),rho(:,3));
    
    [mid_vec13, dec13] =  perp_vec(rho(:,1),rho(:,10),rho(:,3));
    [mid_vec14, dec14] =  perp_vec(rho(:,1),rho(:,3),rho(:,2));
    [mid_vec15, dec15] =  perp_vec(rho(:,1),rho(:,2),rho(:,4));
    [mid_vec16, dec16] =  perp_vec(rho(:,1),rho(:,4),rho(:,9));
    [mid_vec17, dec17] =  perp_vec(rho(:,1),rho(:,9),rho(:,8));
    [mid_vec18, dec18] =  perp_vec(rho(:,1),rho(:,8),rho(:,10));
    [mid_vec19, dec19] =  perp_vec(rho(:,4),rho(:,2),rho(:,7));
    
    
    for i = 1:1
        for j = 1:1
            if (dot(sigma(:,j),mid_vec1)- dec1<= 1e-15) && (dot(sigma(:,j),mid_vec2)- dec2<= 1e-15) && (dot(sigma(:,j),mid_vec3)- dec3<= 1e-15) && (dot(sigma(:,j),mid_vec4)-dec4<= 1e-15) && (dot(sigma(:,j),mid_vec5)- dec5<= 1e-15) && (dot(sigma(:,j),mid_vec6)- dec6<= 1e-15) && (dot(sigma(:,j),mid_vec7)- dec7<= 1e-15)
                %disp('top');
                t=t+1;
            elseif (dot(sigma(:,j),mid_vec8)- dec8<= 1e-15) && (dot(sigma(:,j),mid_vec9)- dec9<= 1e-15) && (dot(sigma(:,j),mid_vec10)- dec10<= 1e-15) && (dot(sigma(:,j),mid_vec11)-dec11<= 1e-15) && (dot(sigma(:,j),mid_vec12)- dec12<= 1e-15)
                %disp('center');
                c=c+1;
            elseif (dot(sigma(:,j),mid_vec13)- dec13<= 1e-15) && (dot(sigma(:,j),mid_vec14)- dec14<= 1e-15) && (dot(sigma(:,j),mid_vec15)- dec15<= 1e-15) && (dot(sigma(:,j),mid_vec16)-dec16<= 1e-15) && (dot(sigma(:,j),mid_vec17)- dec17<= 1e-15) && (dot(sigma(:,j),mid_vec18)- dec18<= 1e-15) && (dot(sigma(:,j),mid_vec19)- dec19<= 1e-15)
               %disp('bottom');
               b=b+1;
            else
                fprintf('%i not done \n', k);
            end
        end
    end
end
t*100/k
c*100/k
b*100/k





function checkpos = is_pos(v,rho,phi)
    res1 = dot(v,rho(:,1))
    res2 = dot(v,rho(:,2))
    res3 = dot(v,rho(:,3))
    res4 = dot(v,rho(:,4))
    res5 = dot(v,phi(:,2))
    res6 = dot(v,phi(:,3))
    if res1 >=0 && res2 >=0 && res3 >=0 && res4 >= 0 && res5 >=0 && res6 >=0
        checkpos = 0;
    elseif (res1 < 0&& res2 < 0 && res6< 0) || (res1 < 0&& res2 < 0 && res3< 0) ||(res1 < 0&& res4 < 0 && res6< 0) || (res1 < 0&& res3 < 0 && res5< 0)|| (res1 < 0&& res4 < 0 && res5< 0)
        checkpos = 1;
    else
        checkpos = 0;
    end
end

function [mid_vec, decider] = perp_vec(r,p,q)
    a = (r(3) - p(3))*(r(2) - q(2)) - (r(2)-p(2))*(r(3)-q(3));
    b = (r(3) - q(3))*(r(1) - p(1)) - (r(3)-p(3))*(r(1)-q(1));
    c = (r(1) - q(1))*(r(2) - p(2)) - (r(1)-p(1))*(r(2)-q(2));    
    d = r(1)*a + r(2)*b + r(3)*c;
    
    mid_vec = [a;b;c];
    decider = d;
end

function pos_ip = pos_ip_with_Xsubset(v)
    if v(1)>=0 && v(2) >=0 && v(3) >=0
        pos_ip = 0;
    elseif v(1)<0 && v(2) >=0 && v(3) >=0 && v(2)+v(3) >= abs(v(1))
        pos_ip = 0;
    else
        pos_ip = 1;
    end
end

%{
function pos_ip = pos_ip_with_Xsubset(v)
    if v(1)>=0 && v(2) >=0 && v(3) >=0
        pos_ip = 0;
    elseif v(1)<0 && v(2) >=0 && v(3) >=0
        if v(2)+v(3) >= abs(v(1))
            pos_ip = 0;
        else
            pos_ip = 1;
        end
    elseif v(1)>=0 && v(2) <0 && v(3) >=0
        %if v(2)==0
        %    pos_ip = 0;
        %else
            pos_ip = 1;
        %end
    elseif v(1)>=0 && v(2) >=0 && v(3) <0
        %if v(3)== 0
        %    pos_ip = 0;
        %else
            pos_ip = 1;
        %end
    elseif v(1)<0 && v(2) <0 && v(3) >=0
        %if v(2)== 0 && v(3) >=v(1)
        %    pos_ip = 0;
        %else
            pos_ip = 1;
        %end
    elseif v(1)<0 && v(2) >=0 && v(3) <0
        %if v(3)== 0 && v(2) >=v(1)
        %    pos_ip = 0;
        %else
            pos_ip = 1;
        %end
    elseif v(1)>=0 && v(2) <0 && v(3) <0
        %if v(2)== 0 && v(3)==0
        %    pos_ip = 0;
        %else
            pos_ip = 1;
        %end
    else% v(1)<0 && v(2)<0 && v(3)<0 
        pos_ip = 1;
    end
end
%}

