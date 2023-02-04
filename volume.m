function RR =volume(pts,A)

% Authors: Teimur Aliev
%          Lebedev Physical Institute of the Russian Academy of Science

a = A(1);
b = A(2);
c = A(3);
tol = min(a,min(b,c))*1e-6; % to exlude points x=y=0!!!!! FIXE
tol2 = 0*tol;

Suplement(1,:) = [2,3,1];
Suplement(2,:) = [3,1,2];
Suplement(3,:) = [1,2,3];
Suplement(4,:) = [2,3,1];
Suplement(5,:) = [3,1,2];
Suplement(6,:) = [1,2,3];


n(1,1:3) = [1,0,0];
n(2,1:3) = [0,1,0];
n(3,1:3) = [0,0,1];
n(4,1:3) = [1,0,0];
n(5,1:3) = [0,1,0];
n(6,1:3) = [0,0,1];
P(1,1:3) = [0,0,0];
P(2,1:3) = [0,0,0];
P(3,1:3) = [0,0,0];
P(4,1:3) = [a,b,c];
P(5,1:3) = [a,b,c];
P(6,1:3) = [a,b,c];

m = 3+sign((pts(:,1)+tol2).*(pts(:,1)-tol2-a))+sign((pts(:,2)+tol2).*((pts(:,2)-tol2)-b))+sign((pts(:,3)+tol2).*(pts(:,3)-tol2-c));
[K vol] = convhulln(pts);

R(1,1:3) = [-1000 -1000 -1000];

SS = find(m==0);

if isempty(SS)
else
    R(2:(length(SS)+1),1:3) = pts(SS,:);
end

ins = [0 0 0 a b c];
ing = [1 2 3 1 2 3];

% rare case one point is outside the cube means full figure is outside the cube
% algorithm is to find planes of one figure which sects points of the snd figure - i.e., distances to the surfaces have different signs
% it is easier to find it for the cube six planes
% index {} - surface of cube, Ki - number of triangle having intersections with given edge

inde(1:6,length(K(:,1))) = 0;
for i = 1:6
    nas{i}  = find(mod(sum(sign(pts(K(:,1),ing(i))-ins(i))+sign(pts(K(:,2),ing(i))-ins(i))+sign(pts(K(:,3),ing(i))-ins(i)),2),3));
    indeS(:,3) = round(1.5*(1+sign((pts(K(:,1),ing(i))-ins(i)).*(pts(K(:,2),ing(i))-ins(i))))); % 0 is both 1st and 2nd vertice lay from different sides of plane    
    indeS(:,2) = round(1+sign((pts(K(:,1),ing(i))-ins(i)).*(pts(K(:,3),ing(i))-ins(i))));
    indeS(:,1) = round(0.5*(1+sign((pts(K(:,2),ing(i))-ins(i)).*(pts(K(:,3),ing(i))-ins(i)))));
    inde(i,:) =  indeS(:,3)+indeS(:,2)+indeS(:,1);% it is the index of one vertices oposed to two other for the triangle crossing plane
    % and 1 is other case then take 1 and 3 vertices
end
for i = 1:1:6
    if isempty (nas{i})
    else
        k = nas{i};
        for j = 1:length(k)
            %first AB section: pts(K(k,1),:)  second points: (1-inde(i,k)).*pts(K(k,2),:)+inde(i,k).*pts(K(k,3),:)
            %snd BC section:
            ha = Suplement(i,3);
            M = Suplement(i,1);
            N = Suplement(i,2);
            x0 = A(M);
            y0 = A(N);
            z0 = floor(i/3.01)*A(ha);
            a0 = inde(i,k(j));
            a1 = mod(a0+1-1,3)+1;
            a2 = mod(a0+2-1,3)+1;
            a0 = mod(a0+0-1,3)+1; % 18.01.2022 FIX
            [I1X{i}(j,:),check1] = plane_line_intersect(n(i,:),P(i,:),pts(K(k(j),a0),:),pts(K(k(j),a1),:));
            [I2X{i}(j,:),check2] = plane_line_intersect(n(i,:),P(i,:),pts(K(k(j),a0),:),pts(K(k(j),a2),:));
            [I1,check1] = plane_line_intersect(n(i,:),P(i,:),pts(K(k(j),a0),:),pts(K(k(j),a1),:));
            [I2,check2] = plane_line_intersect(n(i,:),P(i,:),pts(K(k(j),a0),:),pts(K(k(j),a2),:));
            [I1true I2true] = QUADRAD_AND_LINE([I1(M) I1(N)],[I2(M) I2(N)],x0,y0);
            if (isempty(I1true)||(I1true(1)^2+I1true(2)^2)<tol^2)
            else
                R((end+1),[ha M N]) = [z0 I1true(1) I1true(2)];
            end
            if (isempty(I2true)||(I2true(1)^2+I2true(2)^2)<tol^2)
            else
                R((end+1),[ha M N]) = [z0 I2true(1) I2true(2)];
            end
        end
    end
end

bond = [
    0 0 0
    a b c
    a 0 0
    0 b 0
    0 0 c
    a b 0 
    a 0 c
    0 b c];
[K3 vol3] = convhulln(bond);

in = inhull(bond,pts);
qubepoints = find(in>0); % find vertices of qube that is inside the cube
if isempty(qubepoints)
else
    R(end+1:end+length(qubepoints),:) = bond(qubepoints,:);
end

RR = R(2:end,:);