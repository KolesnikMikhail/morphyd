function [vornbmem,vorvxmem] = symetrizer(pos,paq,dt,A)

%The function symetrizer cuts the boundary Voronoi cells along the domain grains
% Authors: Teimur Aliev and Mikhail Kolesnik
%          Lebedev Physical Institute of the Russian Academy of Science
% e-mail: kolesnik.mikhail@gmail.com

a = A(1);
b = A(2);
c = A(3);

% Reflect points that intersects boundary
if (isempty(paq)||isempty(dt))
    dt = delaunayTriangulation(pos(:,1),pos(:,2),pos(:,3));
    [vornb,vorvx] = voronoiDiagram(dt);

    VB1=find(vornb(:,1) > a);
    VB2=find(vornb(:,2) > b);
    VB3=find(vornb(:,3) > c);
    VB4=find(vornb(:,1) < 0);
    VB5=find(vornb(:,2) < 0);
    VB6=find(vornb(:,3) < 0);
    VornBound=unique([VB1; VB2; VB3; VB4; VB5; VB6]); %Voronoi vertices that are outside the domain
    l=1;
    for j=1:(length(pos(:,1)))
        if isempty(intersect(VornBound,vorvx{j}))
        else
            paq(l) = j;
            l=l+1;
        end
    end
end

paq =[paq,length(pos(:,1))];
cpaq = length(paq);

posin = pos;
mem = length(posin(:,1));
posin(mem+6*cpaq,1:3) = 0;
mis = min(a,min(b,c));
NORM = 1E-6;

mirrors  = repmat(pos(paq,:),6,1);
mirrors2 = zeros([6*cpaq, 3]); 
mirrors2(1:cpaq,1) = -pos(paq,1);
mirrors2((cpaq+1):2*cpaq,2) = -pos(paq,2);
mirrors2((2*cpaq+1):3*cpaq,3) = -pos(paq,3);
mirrors2((3*cpaq+1):4*cpaq,1) = a-pos(paq,1);
mirrors2((4*cpaq+1):5*cpaq,2) = b-pos(paq,2);
mirrors2((5*cpaq+1):6*cpaq,3) = c-pos(paq,3);
randvec2 = NORM*mis*rand(6*cpaq,3);
mirrors = mirrors+2*mirrors2+randvec2;

dt.Points((end+1):(end+6*cpaq),:) = mirrors;

[vornbmem,vorvxmem] = voronoiDiagram(dt);