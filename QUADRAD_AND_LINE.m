function [I1true I2true] = QUADRAD_AND_LINE(I1,I2,x0,y0)

% Authors: Teimur Aliev
%          Lebedev Physical Institute of the Russian Academy of Science

M = 1; N=2;
I1true = [];
I2true = [];

a1 = sign(I1(M)*(x0-I1(M)))+sign(I1(N)*(y0-I1(N))); 
a2 = sign(I2(M)*(x0-I2(M)))+sign(I2(N)*(y0-I2(N)));
if (2-a1>0)&&(2-a2>0) % all lies outside        
    has = I1(M)-I1(N)*(I2(M)-I1(M))/(I2(N)-I1(N));
    if (has*(x0-has)>0)&&(sign((has-I1(M))*(has-I2(M)))<0)
        I1true = [has 0];
    end
    has = I1(M)-(I1(N)-y0)*(I2(M)-I1(M))/(I2(N)-I1(N));
    if (has*(x0-has)>0)&&(sign((has-I1(M))*(has-I2(M)))<0)
        if(isempty(I1true)) 
            I1true = [has y0]; 
        else
            I2true = [has y0];
        end
    end
    has = I1(N)-I1(M)*(I2(N)-I1(N))/(I2(M)-I1(M));
    if (has*(y0-has)>0)&&(sign((has-I1(N))*(has-I2(N)))<0)
        if(isempty(I1true)) 
            I1true = [0 has]; 
        else
            I2true = [0 has];
        end
    end
    has = I1(N)-(I1(M)-x0)*(I2(N)-I1(N))/(I2(M)-I1(M));
    if (has*(y0-has)>0)&&(sign((has-I1(N))*(has-I2(N)))<0)
        if(isempty(I1true)) 
            I1true = [x0 has]; 
        else
            I2true = [x0 has];
        end
    end
else
    if (a1==2)&&(a2==2)
        I1true = I1;
        I2true = I2;
    else
        if (a1==2)
            I1true = I1;
        else
            I1true = I2;
        end
        has = I1(M)-I1(N)*(I2(M)-I1(M))/(I2(N)-I1(N));
        if (has*(x0-has)>0)&&(sign((has-I1(M))*(has-I2(M)))<0)
            I2true = [has 0];
        end
        has = I1(M)-(I1(N)-y0)*(I2(M)-I1(M))/(I2(N)-I1(N));
        if (has*(x0-has)>0)&&(sign((has-I1(M))*(has-I2(M)))<0)
            I2true = [has y0];
        end
        has = I1(N)-I1(M)*(I2(N)-I1(N))/(I2(M)-I1(M));
        if (has*(y0-has)>0)&&(sign((has-I1(N))*(has-I2(N)))<0)
            I2true = [0 has];
        end
        has = I1(N)-(I1(M)-x0)*(I2(N)-I1(N))/(I2(M)-I1(M));
        if (has*(y0-has)>0)&&(sign((has-I1(N))*(has-I2(N)))<0)
            I2true = [x0 has];
        end
    end
    %  R(end+1,[ha M N]) = [z0 P(1) P(2)];
end

if (0)
    figure(1)
    plot([I1(1) I2(1)],[I1(2) I2(2)],'bo-')
    hold on;
    if (isempty(I1true))
    else
        plot([I1true(1) I2true(1)],[I1true(2) I2true(2)],'ro-')
    end
    clear P;
    plot([0 x0 x0 0 0],[0 0 y0 y0 0],'k--')
    axis ([-0.1 1 -0.1 1])
end