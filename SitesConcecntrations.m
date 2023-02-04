function CsI = SitesConcecntrations(pos,Ci,Nvalue,Cs)

%The function SiteConcentrations calculate hydrogen concentration in solid
%solution for each nucleation site, which is equal to the solid solution
%concentration in its Voronoi cell
%
% Authors: Teimur Aliev and Mikhail Kolesnik
%          Lebedev Physical Institute of the Russian Academy of Science
% e-mail: kolesnik.mikhail@gmail.com

global X Y Z n m o correspondance

dx=X/n;
dy=Y/m;
dz=Z/o;

if isempty(Ci)
    CsI=zeros(n*m*o,1)+Cs(end);
else
    for i=1:length(Nvalue)
        if (real(Nvalue(i)) == 0)
            ii = correspondance(i,1);
            jj = correspondance(i,2);
            kk = correspondance(i,3);
            dist=((dx*(ii-0.5)-pos(:,1)).^2+(dy*(jj-0.5)-pos(:,2)).^2+(dz*(kk-0.5)-pos(:,3)).^2).^.5;  %distances from the nucleation site till each hydride
            nearest=find(~(dist-min(dist)));   %the number of the nearest hydride
            CsI(i)=Ci(nearest);
        end
    end
end

CsI=CsI';