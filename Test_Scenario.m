% TEST SCENARIO AND POSTPROCESSING
%
%The scenario assumes cooling of a specimen with a constant cooling rate.
%Initiallyб all hydrogen is in solid solution at temperature above the
%terminal solid solubшlity.
%At the end of the simulation, the most hydrogen is in hydride phase
%
% Authors: Teimur Aliev and Mikhail Kolesnik
%          Lebedev Physical Institute of the Russian Academy of Science
% e-mail: kolesnik.mikhail@gmail.com

clear all;
global X Y Z n m o h ppm correspondance
%_________________________________________________________________________________________________________________________________
%MODEL'S PARAMETERS
dispers = 0.33;    %variance of defect energy distribution
Gdef0 = 8*1.6e-19; %[j] mathematical mean of defect energy distribution
h = 1e-6;          %[m] thickhess of hydrides
ppm = 3.89e24;     %conversion factor from [ppm] to [m^-3]
%_________________________________________________________________________________________________________________________________
%TEMPERATURE SCENARIO
CR = 2/60;       %[C/sec] cooling rate
Tmax = 325+273;   %maximal temperature
Tmin = 320+273;   %minimal temperature
time = 0:((Tmax-Tmin)/CR/1e2):(Tmax-Tmin)/CR; %time steps
Temp = Tmax-time/time(end)*(Tmax-Tmin);   %temperature at time steps
Cs(1) = 200*ppm;  %initial hydrogen consentration in solid solution
sigma = 150e6;     %[Pa] external tension
%_________________________________________________________________________________________________________________________________
%INITIATION OF NUCLEATION SITES AND THE GEOMETRY OF THE DOMAIN (SPECIMEN)
X = .5e-3; %[m] domain length
Y = .5e-3; %[m] domain width
Z = .5e-3; %[m] domain height

n = 50;    %number of nucleation sites per length
m = 50;    %number of nucleation sites per width
o = 50;    %number of nucleation sites per height

dx = X/n;  %[m] nucleation sites grid spacing in length
dy = Y/m;  %[m] nucleation sites grid spacing in length
dz = Z/o;  %[m] nucleation sites grid spacing in length

Nvalue = zeros (n*m*o,1); %matrix of complex values characterizing the hydride structure
%!! Each element of Nvalue corresponds to each nucleation site.
%!! Real part of  Nvalue element defines the hydride radius (positive for radial, negative for tangential, zero for a vacant nucleation site)
%!! Imaginary part of  Nvalue element containes the defect energy (defined below)
correspondance = zeros (n*m*o,3); %contains indexes of 3D matrix of nucleation sites as a column (defined below)
numas = 0;
for i=1:n
    for j=1:m
        for k=1:o
            numas = numas+1;
            correspondance(numas,:) = [i j k];
            Nvalue(numas) = Gdef0*(1+dispers*randn)*sqrt(-1);  %definition of a random defect energy with a given distribution (normal distribution with predefined parameters)
        end
    end
end
%__________________________________________________________________________________________________________________________________
%VARIABLES INITIATION
Active(1) = 0;    %index of activated nucleation sites
pos = [];         %[m] massive of position vectors of activated nucleation sites
paq = [];         %index of Voronoi cells with intersection with boundary  
Ci  = [];         %[m^-3] hydrogen concentrations in solid solution in every Voronoi cell
Vi  = [];         %[m^3] volume of every Voronoi cell
vornbmem = [];    %one of the voronoi results
delaune  = [];    %triangulation
starterpoint = 5; %minimal hydride numbers to start the  Voronoi tessellation 
lentactive = 0;   %last activated center in Active string
%_________________________________________________________________________________________________________________________________
%THE MAIN CONTROL FLOW
for i=1:(size(time,2)-1)  %loop on time
    
    [Cs(i+1),Nvalue,Active,Ci] = NUCLGROWTH((time(i+1)-time(i)),Cs(i),(Temp(i+1)+Temp(i))/2,sigma,Nvalue,Active,Ci,Vi,pos); %nucleation and growth of hydrides
    
    if (isempty(Active))
        AC(i) = 0;
    else %possible nucleation of a new hydride / new hydrides
        if (Active(1)==0)  % cut first artificial null after first activation
            Active(1) = [];
        end
        AC(i) = length(Active); %a number of hydrides
        AC0 = sign(i-1)*AC(max(1,i-1));  % 0 if i=1 and AC(i-1) else  , i.e. the number of hydrides at a previous time step
        if (AC(i) > AC0) %if new hydrides nucleated 
            for IJ = 1:(AC(i)-AC0)  %loop on every nucleated hydride at the current time step
                ii = correspondance(Active(AC0+IJ),1);
                jj = correspondance(Active(AC0+IJ),2);
                kk = correspondance(Active(AC0+IJ),3);
                new = [dx*(ii-0.5) dy*(jj-0.5) dz*(kk-0.5)]+0.01*min(dx,min(dy,dz))*(0.5-rand(1,3));    %new at N(ii,jj,kk);
                if isempty(pos)  %nucleation of the first hydride
                    pos(1,:) = new;
                end
                if (length(pos) < starterpoint) %no Voronoi tessellation (too low number of hydrides)
                    pos(AC0+IJ,1:3) = new; 
                else  %Voronoi tessellation / re-tessellation after adding a new hydride
                    if isempty(vornbmem) %initialize the first tessellation
                        pos(AC0+IJ,1:3) = new; 
                        [vornbmem,vorvxmem] = symetrizer(pos,paq,delaune,[X Y Z]);
                        Ci = Cs(i)+0*(pos(:,1));
                        for xi = 1:length(pos(:,1))
                            [Knull Vi(xi) PONTS] =  V_for_ith(xi,vornbmem,vorvxmem,[X Y Z]); % calculate volumes of Voronoi cells                        
                        end
                    else %re-tessellation after adding a new hydride
                        [pos,Ci,Vi,vornbmem,vorvxmem,dt] = new_point_generator(pos,paq,new,Ci,Vi,vornbmem,vorvxmem,delaune,[X Y Z]); %re-tessellation with a new point
                        for j = 1:1:length(pos(:,1))
                            indu(j) = vorvxmem{j}(1);  
                        end
                        paq = find(indu == 1);
                    end
                end  %end for a Voronoi tessellation / re-tessellation
            end %end for IJ
        end %end for if AC(i)>AC0 - new hydrdies nucleated
    end % if Active is not empty
    fprintf('progress:');disp(i/length(time));
end
%__________________________________________________________________________________________________________________________________
%POSTPROCESSING
%Voronoi tessellation plotting
figure; hold on; title('Voroni tessellation');
for xi = 1:1:length(pos(:,1))
    [Knull Vi(xi) PONTS] =  V_for_ith(xi,vornbmem,vorvxmem,[X Y Z]);  %calculate all volumes
    trisurf(Knull,PONTS(:,1),PONTS(:,2),PONTS(:,3),'FaceColor',rand(1,3),'FaceAlpha',0.5,'EdgeAlpha',1)
end
axis equal;
view([-40 30]);
%__________________________________________________________________________________________________________________________________
%Plotting of the 3D domain with hydrides inside
figure; hold on; title('Hydrodes in the specimen');
plot3([0,X],[0,0],[0,0],'k');  %plotting of the specimen's edges
plot3([0,X],[Y,Y],[0,0],'k'); plot3([0,X],[0,0],[Z,Z],'k'); plot3([0,X],[Y,Y],[Z,Z],'k');  %plotting of the specimen's edges
plot3([0,0],[0,Y],[0,0],'k'); plot3([X,X],[0,Y],[0,0],'k'); plot3([0,0],[0,Y],[Z,Z],'k'); plot3([X,X],[0,Y],[Z,Z],'k');  %plotting of the specimen's edges
plot3([0,0],[0,0],[0,Z],'k'); plot3([X,X],[0,0],[0,Z],'k'); plot3([0,0],[Y,Y],[0,Z],'k'); plot3([X,X],[Y,Y],[0,Z],'k');  %plotting of the specimen's edges
t = 0:pi/100:2*pi;
for i=1:length(pos)
    rhyd=Nvalue(Active(i));
    if rhyd > 0 %plotting of a radial hydride (blue)
        plot3([pos(i,1),pos(i,1)],[pos(i,2)-rhyd,pos(i,2)+rhyd],[pos(i,3),pos(i,3)],'b');
        plot3([pos(i,1),pos(i,1)],[pos(i,2),pos(i,2)],[pos(i,3)-rhyd,pos(i,3)+rhyd],'b');
        plot3(pos(i,1)+t*0,pos(i,2)+rhyd*sin(t),pos(i,3)+rhyd*cos(t),'b');
    else    %plotting of a tangential hydride (red)
        plot3([pos(i,1)-rhyd,pos(i,1)+rhyd],[pos(i,2),pos(i,2)],[pos(i,3),pos(i,3)],'r');
        plot3([pos(i,1),pos(i,1)],[pos(i,2),pos(i,2)],[pos(i,3)-rhyd,pos(i,3)+rhyd],'r');
        plot3(pos(i,1)+rhyd*sin(t),pos(i,2)+t*0,pos(i,3)+rhyd*cos(t),'r');
    end
end
axis equal;
view([-40 30]);
%__________________________________________________________________________________________________________________________________
%plotting of the cross-section at half height
figure; hold on;
title('Cross-section at half height');
xlabel('Tangential direction');
ylabel('Radial direction');
rectangle('Position',[0 0 X Y]);
for i=1:length(pos)
    r_cs=abs(pos(i,3)-round(o/2)*dz);  %the distance from the current hydride center to the cross-section plane
    if r_cs < abs(real(Nvalue(Active(i))))
        rhyd=(real(Nvalue(Active(i)))^2-r_cs^2)^.5;    %hydride radius in the cross-section plane
        if Nvalue(Active(i)) > 0 %plotting of a radial hydride (blue)
            plot([pos(i,1)-h/2,pos(i,1)-h/2],[pos(i,2)-rhyd,pos(i,2)+rhyd],'b');
            plot([pos(i,1)-h/2,pos(i,1)+h/2],[pos(i,2)-rhyd,pos(i,2)-rhyd],'b');
            plot([pos(i,1)-h/2,pos(i,1)+h/2],[pos(i,2)+rhyd,pos(i,2)+rhyd],'b');
            plot([pos(i,1)+h/2,pos(i,1)+h/2],[pos(i,2)-rhyd,pos(i,2)+rhyd],'b');
        else  %plotting of a tangential hydride (red)
            plot([pos(i,1)-rhyd,pos(i,1)-rhyd],[pos(i,2)-h/2,pos(i,2)+h/2],'r');
            plot([pos(i,1)-rhyd,pos(i,1)+rhyd],[pos(i,2)-h/2,pos(i,2)-h/2],'r');
            plot([pos(i,1)-rhyd,pos(i,1)+rhyd],[pos(i,2)+h/2,pos(i,2)+h/2],'r');
            plot([pos(i,1)+rhyd,pos(i,1)+rhyd],[pos(i,2)-h/2,pos(i,2)+h/2],'r');
        end
    end
end
axis equal;
%__________________________________________________________________________________________________________________________________