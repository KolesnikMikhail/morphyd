function [K, V, R] = V_for_ith(is,vornbmem,vorvxmem,A)

% Authors: Teimur Aliev
%          Lebedev Physical Institute of the Russian Academy of Science

tol = min(A(1),min(A(2),A(3)))*1e-6; % to exlude points x=y=0!!!!! FIXE

for i = 1:1:length(is)
    BBC = vorvxmem(is);
    pts = vornbmem(BBC{1}(1:end),:); %? vortexes of k'th neighbour
    in = 3+sign((pts(:,1)+tol).*(pts(:,1)-tol-A(:,1)))+sign((pts(:,2)+tol).*(pts(:,2)-tol-A(:,2)))+sign((pts(:,3)+tol).*(pts(:,3)-tol-A(:,3))); % 
    insideout = 0; % 0 - boundary intersection  1 - no boundary intersection
    if (sum(in)==0) % no intersections
        R = pts;
        insideout = 1;
    else
        if (isinf(pts(1,1)))
            R =  volume(pts(2:end,:),A); % find intersections with cube
        else
            R =  volume(pts(1:end,:),A); % find intersections with cube
        end
    end
    [K V(i)] = convhulln(R,{'QJ','Pp'});       
end