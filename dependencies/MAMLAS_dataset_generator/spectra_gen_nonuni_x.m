% THIS CODE GENERATES ABSORPTION SPECTRA IN SPECIFIED FREQUENCY RANGE
% OF ANY MOLECULE AT ANY THERMODYNAMIC CONDITION OF INTETEST.............. 
% INPUTS: THERMODYNAMIC CONDITIONS, SPECTRAL PARAMETERS, PARTITION FUNCTION
%        *******************NOTES***********************
% Exponent to calculate Pressure Shift in LINE 74 can be 0.0 or 0.96
% for comparing spectra with HITRAN or Spectraplot generated data ....
% .....................................................................
% 'importhitran' function in LINE 35 takes time to read '.par' file of 
% large size. In those cases users are requested to generate an excel file 
% of parameters (with proper indentations) in '.csv' format and uncomment 
% LINES 45 to 53....
%......................................................................
% AUTHOR: DIBYA KANTI GOLUI 

function fit=spectra_gen_nonuni_x(T,P,M,x,L,nu,hitemp_data)

% P=1;                 % Pressure (in Bar)
% M=18.01;             % Molecular Mass of Target Species (in gm/mole)
% x=0.0186;            % Mole fraction of Target Species
% L=10;                % Optical Length (in cm)

% Tmin=1500; Tmax=500; 
% L1=2; L2=L-2*L1;     % Parameters for Non-uniform temp distribution (2T)(Line 59/60)


format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectraplot=readmatrix("H2O HITEMP,x=0.0186,T=1500K,P=1atm,L=10cm,simNum0.csv");
% spectraplot_dat=spectraplot(:,1);
% spectraplot_fit=spectraplot(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nu_min=spectraplot_dat(1);                    % Minimum Wavenumber (in cm-1)
% nu_max=spectraplot_dat(end);                  % Maximum Wavenumber (in cm-1)
% dat=spectraplot_dat;                          % Wavenumber Range (in cm-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nu_min=7184;                                % Minimum Wavenumber (in cm-1)
% nu_max=7185;                                % Maximum Wavenumber (in cm-1)
 dat=nu;             % Wavenumber Range (in cm-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fil1=fopen('q1.txt','r');                   % Partition Function
partfun=fscanf(fil1,'%f %f', [2 Inf]); 
partfun=partfun';
fclose(fil1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lines=importhitran('laser2.par');           % Input in '.par' format taken 
%                                             % directly from HITRAN site
% nu=lines.transitionWavenumber;              % Center-Wavenumber
% s0=lines.lineIntensity;                     % Reference Linestrength
% g_a0=lines.airBroadenedWidth;               % Air-broadened Width
% g_s0=lines.selfBroadenedWidth;              % Self-broadened Width
% Edd=lines.lowerStateEnergy;                 % Lower State Energy
% n=lines.temperatureDependence;              % Temperature Exponent
% d0=lines.pressureShift;                     % Pressure Shift Coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lines=csvread(hitemp_data);                % Input in '.csv' format  
lines=hitemp_data.lines;                    % modified in Excel
nu=lines(:,3);                              % Center-Wavenumber
s0=lines(:,4);                              % Reference Linestrength
g_a0=lines(:,6);                            % Air-broadened Width
g_s0=lines(:,7);                            % Self-broadened Width
Edd=lines(:,8);                             % Lower State Energy
n=lines(:,9);                               % Temperature Exponent
d0=lines(:,10);                             % Pressure Shift Coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Tnon=T.T;
x=x.X;

y=linspace(0,L,length(Tnon))';                        % Discretization of Domain  
dy=L/length(Tnon);
%Tnon=Tmin+(Tmax-Tmin)*(y>L1&y<L1+L2);        % 2-T profile
%Tnon=linspace(Tmin,Tmax,100)';
% figure("Name","first");
% plot(Tnon)
% hold on
% Tnon=Tmax+(Tmin-Tmax)*(y-L/2).^2/(L^2/4);    % Parabolic profile 
S0=s0.*2.479371939e19*P.*x';                   % Unit Conversion of linestrength
                                             % Effects of P,chi are
                                             % considered here

s=ones(length(s0),length(Tnon));
%d=ones(length(s0),length(Tnon));
%g_a=ones(length(s0),length(Tnon));
%g_s=ones(length(s0),length(Tnon));
%nu0=ones(length(s0),length(Tnon));
%delnu_g=ones(length(s0),length(Tnon));
%delnu_l=ones(length(s0),length(Tnon));

d=d0.*(296./Tnon').^0.96;
g_a=g_a0.*(296./Tnon').^n;
g_s=g_s0.*(296./Tnon').^0.75;  % Self-broadened Width 
nu0=nu+P.*(1-x').*d;        % Shifted Center-Wavenumber
s=linestrength(S0,nu0,Edd,Tnon,partfun)*dy;
                                                % Linestrength
delnu_g=0.5*7.162242257e-7.*nu0.*((Tnon'./M).^0.5);
                                                % Gaussian HWHM
delnu_l=P.*(x'.*g_s+(1-x').*g_a);
                                                % Lorentzian HWHM  
%for i=1:length(s0)
    %d(i,:)=d0(i).*(296./Tnon).^0.96;
    %g_a(i,:)=g_a0(i).*(296./Tnon).^n(i);
    %g_s(i,:)=g_s0(i).*(296./Tnon).^0.75;  % Self-broadened Width 
    %nu0(i,:)=nu(i)+P.*(1-x).*d(i,:);        % Shifted Center-Wavenumber
    %s(i,:)=linestrength(S0(i),nu0(i,:),Edd(i),Tnon,partfun)*dy;
                                                % Linestrength
    %delnu_g(i,:)=0.5*7.162242257e-7.*nu0(i,:).*((Tnon'./M).^0.5);
                                                % Gaussian HWHM
    %delnu_l(i,:)=P.*(x.*g_s(i,:)+(1-x).*g_a(i,:));
                                                % Lorentzian HWHM                                            
%     for j=1:length(Tnon)                       
%         d(i,j)=d0(i).*(296/Tnon(j)).^0.96;      % Pressure Shift
%         g_a(i,j)=g_a0(i).*(296/Tnon(j)).^n(i);  % Air-broadened Width
%         g_s(i,j)=g_s0(i).*(296/Tnon(j)).^0.75;  % Self-broadened Width 
%         nu0(i,j)=nu(i)+P.*(1-x).*d(i,j);        % Shifted Center-Wavenumber
%         s(i,j)=linestrength(S0(i),nu0(i,j),Edd(i),Tnon(j),partfun)*dy;
%                                                 %Linestrength
%         delnu_g(i,j)=0.5*7.162242257e-7.*nu0(i,j).*sqrt(Tnon(j)/M);
%                                                 %Gaussian HWHM
%         delnu_l(i,j)=P.*(x.*g_s(i,j)+(1-x).*g_a(i,j));
%                                                 %Lorentzian HWHM
%    end
%end



fit=zeros(length(dat),1);
fit=s;
dat=nu;
% for i=1:length(Tnon)
%     fit=fit+voigt(dat,[nu0(:,i) s(:,i) delnu_g(:,i) delnu_l(:,i)]');                                             
% end


% figure("Name","second");
% plot(dat,fit)  
% title('Absorption Spectrum')
% xlabel('Wavenumber (cm^{-1})')
% ylabel('Absorbance')
% axis tight
% hold on
% legend
% plot(spectraplot_dat,spectraplot_fit,'b')
% hold on
% plot(dat,fit,'r')  
% title('Absorption Spectrum')
% xlabel('Wavenumber (cm^{-1})')
% ylabel('Absorbance')
% axis tight
% hold off
% 
% error=abs(fit-spectraplot_fit)./spectraplot_fit*100;
% plot(dat,error)
% title("Absolute error (%) vs")
% xlabel('Wavenumber (cm^{-1})')
% ylabel('% error')
% axis tight
% disp("mean error: "+string(mean(error))+"  median error: "+ string(median(error))+" max error: "+string(max(error)))
% dat-spectraplot_dat

function q=Q(T,partfun)
%THIS FUNCTION RETURNS VALUE OF PARTITION FUNCTION AT ARBITRARY TEMPERATURE

T=round(T);
T1=T-1;
T2=T+1;
q=zeros(length(T),1);
for iT=1:length(T)
    Tind1=find(partfun(:,1)==T1(iT),1);
    Tind2=find(partfun(:,1)==T2(iT),1);
    q(iT)=interp1(partfun(Tind1:Tind2,1),partfun(Tind1:Tind2,2),T(iT));   
end
% Tind1=find(partfun(:,1)==T1,1);
% Tind2=find(partfun(:,1)==T2,1);
% q=interp1(partfun(Tind1:Tind2,1),partfun(Tind1:Tind2,2),T);
end

function S=linestrength(S0,nu0,Edd,T,partfun)
% THIS FUNCTION RETURNS VALUE OF LINESTRENGTH OF A PARTICULAR TRANSITION
% AT ARBITRARY TEMPERATURE

T0=296;                         % Ref Temperature
h=6.62607004e-34;               % Planck's Constant
c=29979245800;                  % Light Speed at Vacuum
k=1.38064852e-23;               % Boltzmann Constant
S=S0.*((T0./T.*(Q(T0,partfun)./Q(T,partfun)))'.*exp(-h*c.*Edd./k.*(1./T-1/T0)')...
    .*((1-exp(-h*c.*nu0'./(k.*T))))'./(1-exp(-h*c.*nu0'./(k*T0)))');
% S=S0.*(T0./T.*(Q(T0,partfun)/Q(T,partfun)).*exp(-h*c.*Edd/k*(1./T-1/T0))...
%     .*(1-exp(-h*c.*nu0./(k.*T)))./(1-exp(-h*c.*nu0./(k*T0))));
end
end
