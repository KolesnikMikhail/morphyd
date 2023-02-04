function [Cs,Nvalue,activ,Ci] = NUCLGROWTH(time_step,OverSaturation,Temp,sigma,Nvalue,activ,Ci,Vi,pos)

% The function NUCLGROWTH simulates nucleation and growth of hydrides at
% the time step during cooling of an hydrogenated zirconium specimen
%
% Authors: Teimur Aliev and Mikhail Kolesnik
%          Lebedev Physical Institute of the Russian Academy of Science
% e-mail: kolesnik.mikhail@gmail.com

global X Y Z h ppm
tau=time_step;     %[s] time step
Cs=OverSaturation; %[m^-3] mean hydrogen concentration in solid solution
T=Temp;            %[K] temperature
%_________________________________________________________________________________________________________________________________
%PHYSICAL CONSTANTS AND PROEPRTIES OF MATERIALS
R = 8.314463;      %[J/K/mol] Gas constant
kb = 1.380649e-23; %[J/K] Boltzman constant

gamma = 0.065;              %[J/m^2] metal/hydride interface surface energy according to [doi.org/10.1063/1.51021766] and [doi.org/10.1016/j.commatsci.2009.05.025]
vH = 1.7e-29;               %[m^3] atomic volume of hydrogen in hydride according to [doi.org/10.1007/978-1-4471-4195-2]
epsilon = 0.1635;           %volume dilatation in Zr -> ZrH phase transition according to [doi.org/10.1007/978-1-4471-4195-2]
Xa = 1.6/2.6;               %molar fraction of hydrogen in hydride (ZrH1.6)

c1 = 66440; c2=-29630;      %parameters for the solubility precipitation TSSP of type c1*exp(c2/R/T) [ppm] in Zricaloy-4 according to [doi.org/10.1016/j.jnucmat.2011.11.009]
TSS = c1*exp(c2/R./T)*ppm;  %[m^-3] termianl solid solubility (for precipitation) of hydrogen in zirconium

D0 = 7e-7; Q=5360;          %parameters for the hydrogen diffusivity correlation of type D0*exp(-Q/T) [m^2/s] according to [doi.org/10.1016/0022-3115(72)90065-7]
DH = D0*exp(-Q/T);          %[m^2/s] hydrogen diffusion coefficient in zirconium

f0=5000;                    %correlation parameter for the function Fr (radial hydride fraction) - any alloy
Omega = 4.1e-28; ftex = 0.3;%correlation parameter for the function Fr (radial hydride fraction) - Zircaloy-4 (based on data [doi.org/10.1007/BF02652294])
%Omega = 4.5e-28; ftex = 0; %correlation parameter for the function Fr (radial hydride fraction) - Zr-2.5Nb (based on data [doi.org/10.1016/0022-3115(75)90132-4])
Fr=(1-ftex)./(1+f0*exp(-sigma*Omega/kb/T))+ftex;  %radial hydride fraction

Pr = 0.34;                  %Poisson's ratio of zirconium
sy = (590-0.7*T)*1e6;       %[Pa] yield stress of Zircaloy-4 with 0-1400ppm of hydrogen according to [doi.org/10.1016/j.net.2019.07.032]
G = (35.5-0.0198*(T-273))*1e9;  %[Pa] shear modulus of Zircaloy-2 in the range 293 - 773K according to [doi.org/10.1016/0022-3115(75)90071-9]
%_________________________________________________________________________________________________________________________________
%Calculation of cut-off radii for the diffusion task
b = h*epsilon/3;
S = G*vH*b/3/pi*(1+Pr)/(1-Pr);
A = pi*(1+Pr)/4;
TSSP = TSS*exp(-sy*vH*epsilon*A/kb/T);
Dspp = min((X*Y*Z)^(1/3),(nnz(real(Nvalue))/X/Y/Z)^(-1/3)); %[m] mean distance between hydrides
R2 = min(Dspp/2,h/2/vH/(Cs-TSSP)); %[m] outward cut-off radius Rmax
R1 = max(h,S/A/vH/epsilon/sy);     %[m] inward cut-off radius, Rmin
%_________________________________________________________________________________________________________________________________
%THE MAIN CONTROL FLOW
dCs=0;            %[m^-3] change of mean hydrogen concentration in solid solution at the time step
Sink=0;           %[m^-3] total sink of hydrogen from solid solution to hydrides at the time step
if Cs > TSSP 
    if isempty(find(real(Nvalue) ~= 0))
    elseif isempty(Ci)    %before the first Voronoi tessellation
        realNp = find(real(Nvalue)>0);   %radial hydrides growth
        r1 = Nvalue(realNp);
        j1 = 2*pi*DH*(2*pi*r1)*(Cs-TSSP)/(log(R2/R1)+(A*sy*vH*epsilon/kb/T));
        Nvalue(realNp) = Nvalue(realNp)+(j1.*vH/2/pi./r1/h)*tau;  %increment of radial hydrides' radii 
        
        realNm = find(real(Nvalue)<0);  %tangential hydrides growth
        r2 = abs(Nvalue(realNm));
        j2 = 2*pi*DH*(2*pi*r2)*(Cs-TSSP)/(log(R2/R1)+(A*sy*vH*epsilon/kb/T));
        Nvalue(realNm) = Nvalue(realNm)-(j2.*vH/2/pi./r2/h)*tau;  %increment of tangential hydrides' radii
        
        Sink=(sum(j1)+sum(j2))*tau;        
    else                 %growth of hydrides with account of Voronoi tessellation
        
        %growth of radial hydrides
        realNp = activ(ismember(activ, find(real(Nvalue)>0)));
        r1 = Nvalue(realNp);
        j1 = 2*pi*DH*(2*pi*r1).*max(0,Ci(find(ismember(activ, realNp)))-TSSP)./(log(R2/R1)+(A*sy*vH*epsilon/kb/T));
        Ci(find(ismember(activ, realNp)))=Ci(find(ismember(activ, realNp)))-j1*tau./Vi(find(ismember(activ, realNp)))';
        Nvalue(realNp) = Nvalue(realNp)+(j1.*vH/2/pi./r1/h)*tau;  %increment of radial hydrides' radii
                        
        %growth of tangential hydrides
        realNm = activ(ismember(activ, find(real(Nvalue)<0)));
        r2 = abs(Nvalue(realNm));
        j2 = 2*pi*DH*(2*pi*r2).*max(0,Ci(find(ismember(activ, realNm)))-TSSP)./(log(R2/R1)+(A*sy*vH*epsilon/kb/T));
        Ci(find(ismember(activ, realNm)))=Ci(find(ismember(activ, realNm)))-j2*tau./Vi(find(ismember(activ, realNm)))';
        Nvalue(realNm) = Nvalue(realNm)-(j2.*vH/2/pi./r2/h)*tau;  %increment of tangential hydrides' radii
        
        Sink=(sum(j1)+sum(j2))*tau;
    end
    
    %nucleation of new hydrides  
    if isempty(Ci)   %before the first Voronoi tessellation, the hydrogen concentration in solid solution for each nucleation site is equal to the mean solid solution concentration 
        CsIA = Cs+zeros(length(Nvalue),1);
        ActSites = find((real(Nvalue)==0));
        CsI = CsIA(ActSites);
    else             %after the Voronoi tessellation, the hydrogen concentration in solid solution for a nucleation site equal to the H concentration in its Voronoi cell
        CsIA = SitesConcecntrations(pos,Ci,Nvalue,Cs);  %definition of a concentration in solid solution for each nucleation site
        ActSites1 = find((real(Nvalue)==0));
        ActSites2 = find((CsIA > TSSP));
        ActSites = intersect(ActSites1, ActSites2');   %to ignore nucleation sites where concentration of H in solid solution is below TSSP
        CsI = CsIA(ActSites);
    end
    
    Gdef = max(0,imag(Nvalue(ActSites)));
    rn = 2*gamma./(kb*T/vH*Xa*log(CsI/TSSP)-sy*epsilon);  %[m] radius of critical nucleus (for each nucleation site)
    Vn = 4*pi/3*rn.^3;                                    %[m^3] volume of a critical nucleus
    Gibbs = -Vn/vH*Xa*kb*T.*log(CsI/TSSP)+Vn*sy*epsilon+4*pi*rn.^2*gamma;  %[J] Gibbs energy change for critical nucleus
    Zel = (Gibbs/3/pi/kb/T./(Vn/vH).^2).^.5;              %Zeldovich factor, Eq.(2.107) in [doi.org/10.1201/b15829]
    jn = 2*pi*DH*rn.*(CsI-TSSP);                          %[s^-1] diffusion flux to a critical nucleus           
    Rate = Zel.*jn.*exp(-Gibbs/kb/T).*exp(Gdef/kb/T);     %[m^-3*s^-1] nucleation frequency (for each nucleation site)
    
    rand1 = rand(length(ActSites),1);
    rand1c = find(rand1 < Rate*tau);    %random nucleation at each nucleation site
    if ~isempty(rand1c)
        activ(end+1:end+length(rand1c)) = ActSites(rand1c);
        Nvalue(ActSites(rand1c)) = rn(rand1c).*sign(rand(length(rand1c),1)-(1-Fr));
        dCs = pi*sum(rn(rand1c).^2)*h/vH/X/Y/Z;
    end
end
Cs=Cs-dCs-Sink/X/Y/Z;