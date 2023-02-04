function [pos2,Ci2,Vi2,vornbmem2,vorvxmem2,dt] = new_point_generator(pos,paq,new, Ci,Vi,vornbmem,vorvxmem,dt,A)

%The function new_point_generator generates a new Voronoi tessellation after adding a new poit
% Authors: Teimur Aliev and Mikhail Kolesnik
%          Lebedev Physical Institute of the Russian Academy of Science
% e-mail: kolesnik.mikhail@gmail.com

a = A(1); 
b = A(2);
c = A(3);

pos2 = pos;
Ci2 = Ci;
Vi2=Vi;
pos2(end+1,:) = new(:); 
LE = length(pos2(:,1));

if isempty(dt)
    dt = delaunayTriangulation(pos(:,1),pos(:,2),pos(:,3));    
end
dt.Points(end+1,:) = new;
[vornbmem2,vorvxmem2] = symetrizer(pos2,paq,dt,A);

mm = any(dt.ConnectivityList==LE,2);  %non-zero elements of the triangulation which contains the new point
vornc = unique(reshape(dt.ConnectivityList(mm,:),[],1)); % list of neighbour cell of the new point, including itself

for i=1:length(vornc)
    [~, Vi2(vornc(i)), ~] = V_for_ith(vornc(i),vornbmem2,vorvxmem2,A);
end

vornc = vornc(vornc~=LE);  %list of of neighbour cell except the new one

[~, Vi20(vornc), ~] = V_for_ith(vornc,vornbmem2,vorvxmem2,A); %list of neighbour volumes except the new one
[~, Vi(vornc), ~] = V_for_ith(vornc,vornbmem,vorvxmem,A);

Volbefore = sum(Vi(vornc)); 
Volafter = sum(Vi20(vornc)); 
Cbefore = sum(Vi(vornc).*Ci(vornc)');
Cafter = sum(Vi20(vornc).*Ci(vornc)');

Ci2(LE) = (Cbefore-Cafter)/(Volbefore-Volafter);