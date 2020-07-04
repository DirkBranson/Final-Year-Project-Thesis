close all
clear all
clc
%Analytical Design Tool%
%%%%%%%%%%%%%%%% Alex Dunnett - Spherical Revolution %%%%%%%%%%%%%
%Assumptions
%Physical Constants
MaximumCirclePackingDensity = 0.9019; %Max geomtrical packing density for circles, or circle cross-sections.
%-Permativitty
u0 = 4*pi*10^-7; % Pemativitty of Free Space
%-Material Constants
%--Copper
rho_c  = 1.71 *10^-8;%Copper Resistiivty 0.0171 Ohm
rho_C = 8940; %Copper Density kg.m^-3
%--Steel
rho_s = 1.18*10^-7 ; %Steel Resistiviy Ohm.m 
rho_S = 7750; %Steel Density kg.m^-3,
us = 1000; %Relative Permativitty of Steel
%--Aluminium
rho_a = 2.65*10^-8; %Aluminium Resistivity Ohm·m
rho_A = 2700; %kg.m^-3,
%--ABS Plastic
rho_abs = 1.5 * 10^16; %Ohm·m
rho_ABS = 1040; %kg.m^-3,
%--Air
rho_air = 1.5 * 10^16; %Ohm·m
rho_Air = 1.225; %kg.m^-3,

%Design Parameters
%These could be varied and Optimized but are currently static
Stator_Width = 0.06; %Stator Width
Stator_R_In = 0.05; %Stator Inner Radius
Slot_Height = 0.004; %Slot Dimensions
Slot_Width = 0.006; %Slot Dimensions
Phases = 3; %Phases
Stator_Pitch = (2*pi) / 4; %Stators are segments in the Spherical Rotor, i.e. an arclength
P = 4; % Number of Poles - AC Induction Motor
Airgap = 0.002;  %Airgap - Costa and Branco
B_fr = 0; %Bearing Friction,
R_Dia = 0.1; %Rotor Diameter
Shell_Thickness = 0.03; % Costa and Branco
Frequency = 40; %Hz
I_Load = 4; %Amps, the maximum rated current output of the Arudino Motor Shield Rev 3
%I_Load = 7; %Amps, based on the rated current f 
V_Load = 12; %Volts, the maximum rated voltage inut for the Arduino Motor Shield Rev 3
DC_Internal_Resistance = 1.5; %Ohms
Additional_Resistance = 0; %Ohms Additional Resistance to Prevent Overload
Cube_Magnet_Dimensions = 0.005;
ThetaA0 = 0; %Phase Angle of Phase A
Initial_V = 0; %Initial Absolute Velocity of Rotor
%Permant Magnet Motor
Magnet_Poles = 187; %Derived from Python Script "Homogeonous_Sphere"

%Transform from Design Parameters to  Working Parameters
Slot_Distance = sin(Stator_Pitch)*Stator_R_In; %Distance between Slots
StraightWire_Length = Stator_Width; %Force calculation only uses length of straight wire
Coil_Circumference = 2*StraightWire_Length+(2*(pi*Slot_Distance));
R = (R_Dia/2); %Spherical Rotor Radius

%Parameters used by Ma et. al in their paper:
%"A novel method to calculate the thrust of linear inductionmotor based on instantaneous current value" 
%Equivilated to the existing parameters
bn = (2 * pi * R)/(P * Phases); % slot pitch
Kcoil = (Stator_Pitch/(2 * pi)) * P * Phases * 2; % Number of Coils
Tau = (2 * pi * R)/P; %Pole Pitch
P_pair = 2 * P; % p = Number of Pole Pairs
f = Frequency; %f = Frequency
w = f; %angular frequency
Synchronous_V = 2*f*Tau; %Vs = Wave Velocity (tangential)
Relative_V = (Synchronous_V - Initial_V); %Vr_theta = Relative velocity of wave to rotor
Slip = Relative_V / Synchronous_V; %S = Slip
PhaseA = 0; %Relative phase offset of Phase A
t = 0; %Time is presently zero
I_RMS = I_Load/(2^0.5); %RMS current 
b = R + Airgap + (0.5 * Slot_Height); %Centre of rotor, offset from polar origin. Necessary for integrals in Polar Co-ordinates
thetalimit_shell = sin(R/b); % The angle between theta = 0 and the tangent of the rotor shell relative to the origin. Necessary for integrals in Polar Co-ordinates
thetalimit_core = sin((R - Shell_Thickness)/b); % The angle between theta = 0 and the tangent of the rotor core relative to the origin. Necessary for integrals in Polar Co-ordinates
Theta_Z = 0; % Used to distinguish between tangential and radial flux
Theta_Y = pi/2; % Used to distinguish between tangential and radial flux   
MaximumCoilCrossSection_Area = Slot_Height * Slot_Width ; %Stator Slot Area
% Key Design Parameters
Scale = linspace(0.1,10,100); %100 intervals of scale factor. Resulting in rotor diameter from 10mm to 1000mm.
i = 1; %Core Material
j = 1; %Shell Material
AWG = [11 12 13 14 15 16 17 18 19 20 21 22 23];
AWG_WireDiameters = [0.00230 0.00205 0.00182 0.00162 0.00145 0.00129 0.00115 0.001024 0.000912 0.000812 0.000723 0.000644 0.000574]; % m, RS Pro AWG 16
AWG_RatedCurret = [30 25 22.5 20 19 18 16.5 15 14 12.5 11 9 7 4];
%Enamelled Copper Wire is found with American Wire Gauge diameters.
Core_Composition = {'Copper', 'Steel', 'Aluminium', 'ABS','Air'}; %Core Material
Shell_Composition = {'Copper', 'Steel', 'Aluminium', 'ABS','Air'}; %Shell Material
Density = {rho_C, rho_S,rho_A,rho_ABS,rho_Air}; %Material Densities
Resistivity = {rho_c, rho_s, rho_a, rho_abs,rho_air}; %Material Resistivities

AWG_num = size(AWG_WireDiameters,2);
SF_num = size(Scale,2);
Material_Number = size(Core_Composition,2); % Number of Materials
Core_num = Material_Number;
Shell_num = Material_Number;

%Enable I_Load to base current of AWG current operating limits (60degrees)
%Engineering toolbox
for AWG = 1:AWG_num
    WireCrossSectional_Area{AWG} = pi*(AWG_WireDiameters(1,AWG)^2)/4;
    %I_Load{AWG} = AWG_RatedCurret;
end 
%Variation of Wire Gauge
for k = 1:SF_num
   %Variation of Scale Factor, SF denotes array.
    Scale_Factor = Scale(1,k);
    b_SF{k} = b * Scale_Factor;
    Slot_Distance_SF{k} = Slot_Distance * Scale_Factor;
    R_SF{k} = R * Scale_Factor;
    Shell_Thickness_SF{k} = Shell_Thickness * Scale_Factor;
    StatorWidth_SF{k} = 0.06 * Scale_Factor;
    Airgap_SF{k} = Airgap * Scale_Factor;
    StraightWire_Length_SF{k} = StatorWidth_SF{k};
    Coil_Circumference_SF{k} = Coil_Circumference * Scale_Factor;
    Rotor_Diameter_SF{k} = R_SF{k}; 
    MaximumCoilCrossSection_Area_SF{k} = MaximumCoilCrossSection_Area * (Scale_Factor^2);
   
    
    for AWG = 1:AWG_num       
            N_SF{k}{AWG} = round(((MaximumCoilCrossSection_Area_SF{k}/WireCrossSectional_Area{AWG})*MaximumCirclePackingDensity)-0.5); % Number of Turns per Coil
    %Number of coil turns that will fit into stator cross sectional area
    for K = 1:Kcoil
%Slots contribute to the magnetic field by virtue of position as well as
%time.

%Current, Voltage and Resistance of Coils
%Current as a function of position and time (Itk)
            Coil_Resistance_SF{k}{AWG} = (N_SF{k}{AWG} * Coil_Circumference_SF{k} * rho_c) / WireCrossSectional_Area{AWG}; %should probably include Inductance/ Reactance tbh
            XR{k}{AWG} = Coil_Resistance_SF{k}{AWG} + DC_Internal_Resistance + Additional_Resistance;
            Coil_Current_SF{k}{AWG}{K} = V_Load / XR{k}{AWG};
            I_RMS_SF{k}{AWG}{K} = Coil_Current_SF{k}{AWG}{K}/(2^0.5);
            Itk_SF{k}{AWG}{K}=(2^0.5)*N_SF{k}{AWG}*I_RMS_SF{k}{AWG}{K}*sin((-w * t +ThetaA0)+PhaseA+(((2*K)+1)*((bn*pi)/(P_pair*Tau))));%Current in terms of t and k, currently the same everywhere, should be changed per phase
    

    %Approximate indication of Magnetic Field produced by stator windings.
    B_Coil{k}{AWG}{K} = (u0 * Itk_SF{k}{AWG}{K} * us) / (2 * pi); 
    
    %Inductance,Reactance and Impedance
    %L{k}{AWG}{K} = (N_SF{k}{AWG} * B_Coil{k}{AWG}{K})/Itk_SF{k}{AWG}{K};
    L{k}{AWG}{K} = (N_SF{k}{AWG}^2 * WireCrossSectional_Area{AWG} * u0 * us)/StatorWidth_SF{k};
    XL{k}{AWG}{K} = 2 * pi * f * L{k}{AWG}{K};
    Z{k}{AWG}{K} = complex(XR{k}{AWG}, XL{k}{AWG}{K});
    %Preparing variables for iteration. Used to converge on impedance.
    Coil_Current_SF_t{k}{AWG}{K}{1} = Coil_Current_SF{k}{AWG}{K};
    I_RMS_SF_t{k}{AWG}{1} =  I_RMS_SF{k}{AWG};
    Itk_SF_t{k}{AWG}{K}{1} = Itk_SF{k}{AWG}{K};
    B_Coil_t{k}{AWG}{K}{1} = B_Coil{k}{AWG}{K};
    L_t{k}{AWG}{K}{1} = L{k}{AWG}{K} ;
    XL_t{k}{AWG}{K}{1} = XL{k}{AWG}{K};
    Z_t{k}{AWG}{K}{1} = Z{k}{AWG}{K};
    
%         Self - contained iterative process to find impedance
            for iter = 2:30
                Coil_Current_SF_t{k}{AWG}{K}{iter} = V_Load / abs(Z_t{k}{AWG}{K}{iter-1});
                I_RMS_SF_t{k}{AWG}{K}{iter} = Coil_Current_SF_t{k}{AWG}{K}{iter}/(2^0.5);
                Itk_SF_t{k}{AWG}{K}{iter}=(2^0.5)*N_SF{k}{AWG}*I_RMS_SF_t{k}{AWG}{K}{iter}*sin((-w*t+ThetaA0)+PhaseA+(((2*K)+1)*((bn*pi)/(P_pair*Tau))));
                B_Coil_t{k}{AWG}{K}{iter} = (u0 * Itk_SF_t{k}{AWG}{K}{iter} * us) / (2 * pi);
                L_t{k}{AWG}{K}{iter} = (N_SF{k}{AWG} * B_Coil_t{k}{AWG}{K}{iter-1})/Itk_SF_t{k}{AWG}{K}{iter};
                XL_t{k}{AWG}{K}{iter} = 2 * pi * f * L_t{k}{AWG}{K}{iter};
                Z_t{k}{AWG}{K}{iter} = complex(XR{k}{AWG}, XL_t{k}{AWG}{K}{iter});
                Zabs_t{k}{AWG}{K}{iter} = abs(Z_t{k}{AWG}{K}{iter});
            end
        
        Coil_Current_SF{k}{AWG}{K} = Coil_Current_SF_t{k}{AWG}{K}{iter};
        I_RMS_SF{k}{AWG}{K} = I_RMS_SF_t{k}{AWG}{K}{iter};
        Itk_SF{k}{AWG}{K} = Itk_SF_t{k}{AWG}{K}{iter};
        B_Coil{k}{AWG}{K} = B_Coil_t{k}{AWG}{K}{iter};
        L{k}{AWG}{K} = L_t{k}{AWG}{K}{iter};
        XL{k}{AWG}{K} = XL_t{k}{AWG}{K}{iter};
        Z{k}{AWG}{K} = Z_t{k}{AWG}{K}{iter};
        Zabs{k}{AWG}{K} = abs(Z_t{k}{AWG}{K}{iter});
    
        %There are two current limits. One dictated by the motor controller
        %(4A) and one dictated by the breakdown voltage of the enamelled
        %wire being used. The value of I_load reflects this.
        Required_Resistance{k}{AWG}{K} = V_Load / I_Load;
        Resistor{k}{AWG}{K} = Required_Resistance{k}{AWG}{K} - Zabs{k}{AWG}{K};
        I_load_final{k}{AWG}{K} = V_Load./(Zabs{k}{AWG}{K} + Resistor{k}{AWG}{K});
        %Itk_SF{k}{AWG}{K} = I_load_final{k}{AWG}{K}; %current limited to 4A
    end 
     
%All current contribution to flux summed across one stator. Will later  be 
%needed for summing across all stators
    B_Stator{k}{AWG} = 0;
    Itk_Stator_SF{k}{AWG} = 0;
    
        for K = 1:Kcoil
            B_Stator{k}{AWG} = B_Stator{k}{AWG} + B_Coil{k}{AWG}{K};
            Itk_Stator_SF{k}{AWG} = Itk_Stator_SF{k}{AWG} + Itk_SF{k}{AWG}{K};
        end
    end
    for i = 1 : Material_Number
    %Inertia calculation
        Core_Volume_Sphere{k}{i} = (4*pi*(R_SF{k}^3) / 3);
        Shell_Volume_Sphere{k}{i} = (4*pi*((R_SF{k}+Shell_Thickness_SF{k})^3) / 3) - (4*pi*(R_SF{k}^3) / 3);
        Core_Mass_Sphere{k}{i} = Core_Volume_Sphere{k}{i} * Density{i};
        Shell_Mass_Sphere{k}{i} = Shell_Volume_Sphere{k}{i} * Density{i};
        Core_Inertia_Sphere{k}{i} = (2/5) * Core_Mass_Sphere{k}{i} * R_SF{k}^2;
        Shell_Inertia_Sphere{k}{i} = (2/5) * ((Core_Mass_Sphere{k}{i} * (R_SF{k} + Shell_Thickness_SF{k})^2) - (Core_Mass_Sphere{k}{i} * R_SF{k}^2));
        Core_Volume_Cylinder{k}{i} = pi * (R_SF{k}^2) * 2 * (StatorWidth_SF{k}/2) ;
        Shell_Volume_Cylinder{k}{i} = ((pi*(R_SF{k} + Shell_Thickness_SF{k})^2) * 2 * (StatorWidth_SF{k}/2)) - ((pi*(R_SF{k})^2) * 2 * (StatorWidth_SF{k}/2));
        Core_Mass_Cylinder{k}{i} = Core_Volume_Cylinder{k}{i}*Density{i};
        Shell_Mass_Cylinder{k}{i} = Shell_Volume_Cylinder{k}{i}*Density{i};
        Core_Inertia_Cylinder{k}{i} = (1/2) * Core_Mass_Cylinder{k}{i} * R_SF{k}^2;
        Shell_Inertia_Cylinder{k}{i} = ((1/2) * Core_Mass_Cylinder{k}{i} * (R_SF{k}+Shell_Thickness_SF{k})^2) - ((1/2) * Core_Mass_Cylinder{k}{i} * R_SF{k}^2);

            for AWG = 1:AWG_num
                %EMF generated by changing magnetic field
                E_ind{k}{AWG}{i} = B_Stator{k}{AWG} * cos(2 * Theta_Z) * Relative_V * 2 * Stator_Width;
                %Eddy currents assumed to flow along the contours of the electric
                %field flux. No laminations in rotor so there would be significant
                %eddy current losses, but that is assumed not to be the case here.
                dI_ind{k}{AWG}{i} = E_ind{k}{AWG}{i} / ((N_SF{k}{AWG} * Coil_Circumference_SF{k} * Resistivity{i}) / (WireCrossSectional_Area{AWG}));  


                %differential form of force acting over an element of the rotor
                Lorentz_Forces{k}{AWG}{i} = B_Stator{k}{AWG} * dI_ind{k}{AWG}{i} * Stator_Width;

                %Function for integration. A polar co-ordinate system is assumed
                %such that a 2D integration of flux penetration can occur.
                %Direction of flux is approximate by the sin or cos (2.theta0
                %function at the. Radial or Tangential Penetration.
                %A double integral over both theta and radius was necessary.
                %First integral was performed analytical.Second was done with
                %numerical methods.
                Flux_Normalized_Function_Radial_Core{k} = @(theta_variable) (log( b_SF{k} ) - log( b_SF{k}.^2 - 2.*b_SF{k}.*(R_SF{k} - Shell_Thickness_SF{k}).*cos(theta_variable) + (R_SF{k} - Shell_Thickness_SF{k})^2)) .* cos(1 .* theta_variable) ;
                Flux_Normalized_Function_Tangential_Core{k} = @(theta_variable) (log( b_SF{k} ) - log( b_SF{k}.^2 - 2.*b_SF{k}.*(R_SF{k} - Shell_Thickness_SF{k}).*cos(theta_variable) + (R_SF{k} - Shell_Thickness_SF{k})^2)) .* sin(1 .* theta_variable) ;
                Flux_Normalized_Function_Radial_Shell{k} = @(theta_variable) (log( b_SF{k} ) - log( b_SF{k}.^2 - 2.*b_SF{k}.*R_SF{k}.*cos(theta_variable) + R_SF{k}^2)) .* cos(1 .* theta_variable);
                Flux_Normalized_Function_Tangential_Shell{k} = @(theta_variable) (log( b_SF{k} ) - log( b_SF{k}.^2 - 2.*b_SF{k}.*R_SF{k}.*cos(theta_variable) + R_SF{k}^2)) .* sin(1 .* theta_variable);
                %Simple numerical integration
                Flux_Integral_Radial_Core{k}           =integral(Flux_Normalized_Function_Radial_Core{k},-thetalimit_core,thetalimit_core);
                Flux_Integral_Tangential_Core{k}       =integral(Flux_Normalized_Function_Tangential_Core{k},-thetalimit_core,thetalimit_core);
                Flux_Integral_Radial_Shell{k}          =integral(Flux_Normalized_Function_Radial_Shell{k},-thetalimit_shell,thetalimit_shell) - Flux_Integral_Radial_Core{k};
                Flux_Integral_Tangential_Shell{k}      =integral(Flux_Normalized_Function_Tangential_Shell{k},-thetalimit_shell,thetalimit_shell) - Flux_Integral_Tangential_Core{k};
                %Total force acting on Shell and Core of Rotor. Unsurprisingly
                %Tangential contribution is negligible.
                Force_Shell_TANG{k}{AWG}{i} = Lorentz_Forces{k}{AWG}{i} * Flux_Integral_Radial_Shell{k};
                Force_Shell_RADI{k}{AWG}{i} = Lorentz_Forces{k}{AWG}{i} * Flux_Integral_Tangential_Shell{k};
                Force_Core_TANG{k}{AWG}{i} = Lorentz_Forces{k}{AWG}{i} * Flux_Integral_Radial_Core{k};
                Force_Core_RADI{k}{AWG}{i} = Lorentz_Forces{k}{AWG}{i} * Flux_Integral_Tangential_Core{k};

                %Contribution from both tangential and radial is summed.
                Total_Force_Shell_RADI{k}{AWG}{i} = Force_Shell_RADI{k}{AWG}{i}*2*StraightWire_Length_SF{k};
                Total_Force_Shell_TANG{k}{AWG}{i} = Force_Shell_TANG{k}{AWG}{i}*2*StraightWire_Length_SF{k};
                Total_Force_Core_RADI{k}{AWG}{i} = Force_Core_RADI{k}{AWG}{i}*2*StraightWire_Length_SF{k};
                Total_Force_Core_TANG{k}{AWG}{i} = Force_Core_TANG{k}{AWG}{i}*2*StraightWire_Length_SF{k};
                %Performance Characteristics are calculated
                Torque_Shell{k}{AWG}{i} = R_SF{k} * Total_Force_Shell_TANG{k}{AWG}{i};
                Torque_Core{k}{AWG}{i} = (R_SF{k} - Shell_Thickness_SF{k}) * Total_Force_Core_TANG{k}{AWG}{i};
                %Forces is multiplied by number of stators or total stator
                Force_Total{k}{AWG}{i} = ((2*pi)/Stator_Pitch) * (Force_Core_TANG{k}{AWG}{i} + Force_Shell_TANG{k}{AWG}{i});
                Torque_Total{k}{AWG}{i} =((2*pi)/Stator_Pitch) * (Torque_Shell{k}{AWG}{i} + Torque_Core{k}{AWG}{i});
                Inertia_Total_Sphere{k}{i} = Core_Inertia_Sphere{k}{i} + Shell_Inertia_Sphere{k}{i};
                Acceleration_Sphere{k}{AWG}{i} = Torque_Total{k}{AWG}{i}/Inertia_Total_Sphere{k}{i};
            end
     end

end
%Goodness Factor Calc
for k = 1:SF_num
for Shell = 1: Shell_num
Goodness_Factor{Shell}{k} = (w * u0 * Slot_Distance_SF{k}^2)/(Resistivity{Shell} * Airgap_SF{k});
end
end
%Preparing core and shell combinations
        for k = 1:SF_num  
                 for Core = 1:Core_num
                      for Shell = 1:Shell_num
                          Inertia_PTWRK{k}{Core}{Shell} = Core_Inertia_Sphere{k}{Core} + Shell_Inertia_Sphere{k}{Shell};
                            for AWG = 1:AWG_num
                                  
                                   Torque_PTWRK{k}{AWG}{Core}{Shell} = Torque_Core{k}{AWG}{Core} + Torque_Shell{k}{AWG}{Shell};
                                   Force_PTWRK{k}{AWG}{Core}{Shell} = Total_Force_Core_TANG{k}{AWG}{Core} + Total_Force_Shell_TANG{k}{AWG}{Core};
                                   Acceleration_PTWRK{k}{AWG}{Core}{Shell} = Torque_PTWRK{k}{AWG}{Core}{Shell}/Inertia_PTWRK{k}{Core}{Shell};
                            end
                      end
                 end
        end    
        %Creating Scale_Factor Vectors for use in graphs and figures
       for Core = 1:Core_num
           for Shell = 1:Shell_num
               for k = 1:SF_num
                   Inertia_PTWRK_SizeVector{Core}{Shell}(1,k) = Core_Inertia_Sphere{k}{Core} + Shell_Inertia_Sphere{k}{Shell};
                   for AWG = 1:AWG_num
                       Acceleration_PTWRK_SizeVector{AWG}{Core}{Shell}(1,k) = Acceleration_PTWRK{k}{AWG}{Core}{Shell};
                       Torque_PTWRK_SizeVector{AWG}{Core}{Shell}(1,k) = Torque_Core{k}{AWG}{Core} + Torque_Shell{k}{AWG}{Shell};
                       Force_PTWRK_SizeVector{AWG}{Core}{Shell}(1,k) = Total_Force_Core_TANG{k}{AWG}{Core} + Total_Force_Shell_TANG{k}{AWG}{Shell};
                   end
               end
           end
       end
        
        %Additional colours means validation more likely to be successful
          MapColours = colormap(jet(Core_num^2));   
          
        for col = 1:Core_num^2
        Colours{col} = MapColours(col,1:3);
        end
        Colour_Array = reshape(Colours,Core_num,Shell_num);
        
      
      %For size dependant plots
        for k = 1:SF_num
            R_SF_SizeVector(1,k) = R_SF{k};
            Flux_Integral_Radial_Core_SizeVector(1,k) = Flux_Integral_Radial_Core{k};
            Flux_Integral_Tangential_Core_SizeVector(1,k) = Flux_Integral_Tangential_Core{k};
            Flux_Integral_Radial_Shell_SizeVector(1,k) = Flux_Integral_Radial_Shell{k};
            Flux_Integral_Tangential_Shell_SizeVector(1,k) = Flux_Integral_Tangential_Shell{k};
            Lorentz_Forces_SizeVector(1,k) = Lorentz_Forces{k}{AWG}{1};
            dI_ind_SizeVector(1,k) = dI_ind{k}{AWG}{1};
            B_Coil1_SizeVector(1,k) = B_Coil{k}{AWG}{1};
            B_Coil2_SizeVector(1,k) = B_Coil{k}{AWG}{2};
            B_Coil3_SizeVector(1,k) = B_Coil{k}{AWG}{3};
            B_Coil4_SizeVector(1,k) = B_Coil{k}{AWG}{4};
            B_Coil5_SizeVector(1,k) = B_Coil{k}{AWG}{5};
            B_Coil6_SizeVector(1,k) = B_Coil{k}{AWG}{6};
            B_Stator_SizeVector(1,k) = B_Stator{k}{AWG};
            PM_Magnoot(1,k) = 0.2348; %Equivilent field strength at 2mm airgap of N45 Magnets
            E_ind_SizeVector(1,k) = E_ind{k}{AWG}{1};
            for Shell = 1: Shell_num
                Goodness_Factor_SizeVector{Shell}(1,k) = Goodness_Factor{Shell}{k};
            end
        end
        %Sanity Check
        Flux_Integral_Radial_Rotor_SizeVector = Flux_Integral_Radial_Core_SizeVector + Flux_Integral_Radial_Shell_SizeVector;
        
        
        %For size and AWG dependant plots
        %This is important as the number of conductors will be affected.
      for k = 1:SF_num  
        for AWG = 1:AWG_num
               N_SF_Mesh(AWG,k) = N_SF{k}{AWG};
               Itk_Stator_SF_Mesh(AWG,k) = Itk_Stator_SF{k}{AWG};
        end
      end
      
      %For induction iteration and AWG plots
      %Acting over a single coil, in Phase B [K = 2]
for iteration = 1:iter
    for AWG = 1:AWG_num
         Reactance_Iterative(iteration,AWG) = XL_t{10}{AWG}{2}{iteration};
    end
end
      iteration = 1:1:iter;

%Plotting Reactance as a proportion of Impedance for AWG

%Again, Induction calculations take places on coil 2
for AWG = 1:AWG_num
    Zabs_plot(1,AWG) = Zabs{10}{AWG}{2};
    XL_plot(1,AWG) = XL{10}{AWG}{2};
    XR_plot(1,AWG) = XR{10}{AWG};
    XL_over_Zabs(1,AWG) = XL{10}{AWG}{2} / Zabs{10}{AWG}{2};
    XR_over_Zabs(1,AWG) = XR{10}{AWG} / Zabs{10}{AWG}{2};
end




 %Legend creation for all patchwork rotors
       for Core = 1:Core_num
           for Shell = 1:Core_num
               labels{Core}{Shell} = strcat(Core_Composition{Core},' Core &',{' '},Core_Composition{Shell},' Shell');
               label_array{Core,Shell} = labels{Core}{Shell}; 
           end
       end
label_array_reshape = reshape(label_array,1,Core_num^2);
label_array_string = string(label_array_reshape);

           for Shell = 1:Core_num
               labels_shell{Shell} = strcat(Core_Composition{Shell},' Shell');
               label_array_shell{Shell} = labels_shell{Shell}; 
           end
label_array_reshape_shell = reshape(label_array_shell,1,Shell_num);
label_array_string_shell = string(label_array_reshape_shell);

%Figures
  figure('Position', get(0, 'Screensize'))
    grid on
    hold on      
      plot(R_SF_SizeVector, B_Coil1_SizeVector,'r'); 
    plot(R_SF_SizeVector, B_Coil2_SizeVector,'g');
     plot(R_SF_SizeVector, B_Coil3_SizeVector,'b');
      plot(R_SF_SizeVector, B_Coil4_SizeVector,'c');
       plot(R_SF_SizeVector, B_Coil5_SizeVector,'y');
        plot(R_SF_SizeVector, B_Coil6_SizeVector,'m');
%         plot(R_SF_SizeVector, B_Stator_SizeVector,'k');
        plot(R_SF_SizeVector, PM_Magnoot,'k');
      set(gca,'xscale','log')
      xlim([0.01 0.5])
      xlabel('Rotor Radius [m]')
      ylabel('Flux Density [T]')
      legend('Coil 1','Coil 2','Coil 3','Coil 4','Coil 5','Coil 6','Eq. N45 5mm^3 Magnet')
      title('Stator Total and Coil Contribution to Flux Density')
saveas(gcf,'Stator Total and Coil Contribution to Flux Density','png')
hold off  

 figure('Position', get(0, 'Screensize'))
       hold on
       grid
               plot(R_SF_SizeVector,Flux_Integral_Radial_Shell_SizeVector,'c')
               plot(R_SF_SizeVector,Flux_Integral_Radial_Core_SizeVector,'y')
               plot(R_SF_SizeVector,Flux_Integral_Radial_Rotor_SizeVector,'m') 
       title('Magnitude of Flux Normal to Rotor vs Scale Factor')
      ylabel('Dimensionless Normalized Flux Density')
      xlabel('Rotor Radius /m')
      set(gca,'xscale','log')
      xlim([0.005 0.5])
     legend('Shell Radial_Flux','Core Radial_Flux','Rotor Radial Flux');
     saveas(gcf,'Magnitude of Flux Normal to Rotor vs Scale Factor','png')
     hold off
     %AWG = 6 %for graphs
       figure('Position', get(0, 'Screensize'))
       hold on
       grid
       for Core = 1:Core_num
           for Shell = 1:Shell_num
                   plot(R_SF_SizeVector, Torque_PTWRK_SizeVector{AWG}{Core}{Shell},'Color',Colour_Array{Core,Shell})
                   set(gca,'xscale','log')       
           end
       end  
       title('Torque vs Radius for Core and Shell Material Selection')
      ylabel('Torque /Nm')
      xlabel('Rotor Radius /m')
     legend(label_array_string);
saveas(gcf,'Torque vs Radius for Core and Shell Material Selection','png')
hold off

 figure('Position', get(0, 'Screensize'))
       hold on
       grid
       for Core = 1:Core_num
           for Shell = 1:Shell_num
                   plot(R_SF_SizeVector, Acceleration_PTWRK_SizeVector{AWG}{Core}{Shell},'Color',Colour_Array{Core,Shell})
                   set(gca,'yscale','log')       
           end
       end  
       title('Acceleration vs Radius for Core and Shell Material Selection')
      ylabel('Acceleration /ms-2')
      xlabel('Rotor Radius /m')
      set(gca,'xscale','log')
      xlim([0.01 0.5])
      ylim([10^-4 10])
     legend(label_array_string);
saveas(gcf,'Acceleration vs Radius for Core and Shell Material Selection','png')
hold off

 figure('Position', get(0, 'Screensize'))
       hold on
       grid
       for Core = 1:Core_num
           for Shell = 1:Shell_num
                   plot(R_SF_SizeVector, Inertia_PTWRK_SizeVector{Core}{Shell},'Color',Colour_Array{Core,Shell})  
                   set(gca,'xscale','log')
                   set(gca,'yscale','log')
           end
       end  
       title('Inertia vs Radius for Core and Shell Material Selection')
      ylabel('Inertia /kgm-2')
      xlabel('Rotor Radius /m')
     legend(label_array_string);
     xlim([0.01 0.5])
saveas(gcf,'Inertia vs Radius for Core and Shell Material Selection','png')
hold off

 figure('Position', get(0, 'Screensize'))
       hold on
       grid
       for Core = 1:Core_num
           for Shell = 1:Shell_num
                   plot(R_SF_SizeVector, Force_PTWRK_SizeVector{AWG}{Core}{Shell},'Color',Colour_Array{Core,Shell})
                   set(gca,'xscale','log')       
           end
       end  
       title('Force vs Radius for Core and Shell Material Selection')
      ylabel('Force /N')
      xlabel('Rotor Radius /m')
      xlim([0.01 0.5])
     legend(label_array_string);
saveas(gcf,'Force vs Radius for Core and Shell Material Selection','png')
hold off

% Surface plot time
WireCrossSectional_Area = cell2mat(WireCrossSectional_Area);
[Radius, AmericanWireGauge_100] = meshgrid(R_SF_SizeVector,(23+1-AWG_num):23);
[Current,AmericanWireGauge_13] = meshgrid(Itk_Stator_SF_Mesh(:,10),(23+1-AWG_num):23);

      Current(1,:) = Itk_Stator_SF_Mesh(1:AWG,10);
      Turns(1,:) = N_SF_Mesh(1:AWG,10);
      Current_Density(1,:) = Current(1,:)./(WireCrossSectional_Area(1:AWG) .* Turns(1,:))
% Convert cell structures for F. etc. into arrays with indexs (Radius and
% AWG
for AWG = 1:AWG_num
    for Core = 1 : Core_num
        for Shell = 1 : Shell_num
Torque_PTWRK_Surface_Radius{Core}{Shell}(AWG,1:SF_num) = Torque_PTWRK_SizeVector{AWG}{Core}{Shell}(1,:);
Torque_PTWRK_Surface_AWG{Core}{Shell}(AWG,1:AWG_num) = Torque_PTWRK_SizeVector{AWG}{Core}{Shell}(:,10);
        end
    end
end

%Surface variable generator for current density
for k = 1:SF_num
    for AWG = 1 : AWG_num
           Current_Density_Surface(AWG,k) = Itk_Stator_SF{k}{AWG}./(WireCrossSectional_Area(1,AWG) .* N_SF{k}{AWG});
    end
end

figure('Position', get(0, 'Screensize'))
surf(Radius,AmericanWireGauge_100,Current_Density_Surface,'FaceColor','interp','EdgeColor','interp','FaceAlpha',0.5)
xlabel('Radius/m');
ylabel('AWG');
zlabel('Current Density A/m-2'); 
set(gca,'xscale','log')  
saveas(gcf,'Torque vs AWG & Scale Factor - Log','png')
xlim([0.01 0.5])
figure('Position', get(0, 'Screensize'))
surf(Radius,AmericanWireGauge_100,Torque_PTWRK_Surface_Radius{1}{1},'FaceColor','interp','EdgeColor','interp','FaceAlpha',0.5)
xlabel('Radius/m');
ylabel('AWG');
zlabel('Torque /Nm');
set(gca,'xscale','log')    
saveas(gcf,'Torque vs AWG & Scale Factor - Log','png')


       figure('Position', get(0, 'Screensize'))
       plot(AmericanWireGauge_100(:,1),Current_Density(1,:),'r')
       grid on
       title('AWG vs Current')
       xlabel('AWG')
       ylabel('Current Density per Winding /Am-2')
       ylim([0 2*(10^6)])
       saveas(gcf,'Current_Density vs AWG','png')
       
      
       figure('Position', get(0, 'Screensize'))
       plot(AmericanWireGauge_100(:,1),Turns,'r')
       title('AWG vs Turns')
       xlabel('AWG')
       ylabel('Turns per Coil')
       grid on
       saveas(gcf,'AWG vs Turns','png')
       
       
         
           figure('Position', get(0, 'Screensize'))
           hold on
           grid on
           for AWG = 1:AWG_num
          plot(iteration,Reactance_Iterative(:,AWG))
           end
           legend('AWG 11','AWG 12','AWG 13','AWG 14','AWG 15','AWG 16','AWG 17','AWG 18','AWG 19','AWG 20','AWG 21','AWG 22','AWG 23')
          title('Reactance Convergance for all AWG')
          xlabel('Iterations');
          ylabel('Reactance /Ohms')
          saveas(gcf,'Reactance Convergance for all AWG','png') 
          hold off
          
%Inductance and Resistance variation with AWG

figure('Position', get(0, 'Screensize'))
hold on
grid on
plot(AmericanWireGauge_100(:,1),XR_over_Zabs,'r')
plot(AmericanWireGauge_100(:,1),XL_over_Zabs,'g')
xlabel('AWG')
ylabel('Contribution to Impedance')
legend('Resistance as a proportion of Impedance','Reactance as a proportion of Impedance')
title('Proportion of Impedance Contributed by Reactance and Resistance')
hold off
saveas(gcf,'React,Resist,Impede','png')



figure('Position', get(0, 'Screensize'))
hold on
grid on
plot(AmericanWireGauge_100(:,1),Zabs_plot,'r')
xlabel('AWG')
ylabel('Impedance /Ohms')
xlim([11 23])
ylim([0 7])
title('Impedance against AWG')
saveas(gcf,'Impedance against AWG','png')

figure('Position', get(0, 'Screensize'))
hold on
grid on
for Shell = 1 : Shell_num
plot(R_SF_SizeVector,Goodness_Factor_SizeVector{Shell})
end
xlabel('Rotor Radius /m')
ylabel('Goodness Factor')
title('Goodness Factor against Scale Factor')
set(gca,'xscale','log') 
legend(label_array_string_shell)
saveas(gcf,'Goodness Factor against Scale Factor','png')