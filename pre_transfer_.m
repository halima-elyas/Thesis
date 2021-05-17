

clc; clear all; close all;

%--------------------------------------------------------------------------

import org.opensim.modeling.*

% Read in osim model
osimModel=Model('pretransfer.osim');    % create model object by loading the model from file
state=osimModel.initSystem;                      % initialise the model at default state

% Obtain coordinates
coords=osimModel.getCoordinateSet();
flexion=coords.get('flexion');         % flexion
deviation=coords.get('deviation');     % deviation
elbow=coords.get('EL_x');              % elbow flexion 

% Get min/max coordinate values.
rMaxFlex=flexion.getRangeMax();
rMinFLex=flexion.getRangeMin();
rMaxDev=deviation.getRangeMax();
rMinDev=deviation.getRangeMin();
rMaxE=elbow.getRangeMax();
rMinE=elbow.getRangeMin();    

% Get the set of muscles that are in the original model
muscles = osimModel.getMuscles(); 
nMuscles = muscles.getSize();

% Get a handle to the brachioradialis 1, 2 and 3
br1 = osimModel.getMuscles().get("brachiorad_1");
br2 = osimModel.getMuscles().get("brachiorad_2");
br3 = osimModel.getMuscles().get("brachiorad_3");
 
% Get a handle to the ECRB
ECRB = osimModel.getMuscles().get("ECRB");

% Get a handle to rest of wrist extensors
ECU = osimModel.getMuscles().get("ECU");
EDCI = osimModel.getMuscles().get("EDCI");
EDCL = osimModel.getMuscles().get("EDCL");
EDCR = osimModel.getMuscles().get("EDCR");
EDCM = osimModel.getMuscles().get("EDCM");
EDM = osimModel.getMuscles().get("EDM");
EIP = osimModel.getMuscles().get("EIP");
EPL = osimModel.getMuscles().get("EPL");

% Get a handle to rest of elbow flexors
brachialis_1=osimModel.getMuscles().get("brachialis_1");
brachialis_2=osimModel.getMuscles().get("brachialis_2");
brachialis_3=osimModel.getMuscles().get("brachialis_3");
brachialis_4=osimModel.getMuscles().get("brachialis_4");
brachialis_5=osimModel.getMuscles().get("brachialis_5");
brachialis_6=osimModel.getMuscles().get("brachialis_6");
brachialis_7=osimModel.getMuscles().get("brachialis_7");
bic_b_1=osimModel.getMuscles().get("bic_b_1");
bic_b_2=osimModel.getMuscles().get("bic_b_2");

actuator_set=osimModel.getActuators;
force_set=osimModel.getForceSet;
ECRB_force=osimModel.getForceSet.get("ECRB");

%-----------------------------Solution------------------------------------%

% Define solution matrices 
ecrb_f=[];
ecrb_d=[];
ECRB_e=[];
br_e=[];
el_flex=[];
moment_f=[];
moment_d=[];
moment_e=[];
total_m=[];
    for i=0:nMuscles-1
        muscles.get(i).setActivation(state,1.0)
        afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i)); 
    end

    % FLEXION: step through 100 points in ROM: 
    for f=linspace(rMinFLex,rMaxFlex,100)
        flexion.setValue(state, f);

         for i=0:nMuscles-1
            muscles.get(i).setActivation(state,1.0)
            afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
         end

        osimModel.equilibrateMuscles(state); 

        % Compute the moment arm and force at this coordinate.
        % ECRB:
        ECRB_ma=ECRB.computeMomentArm(state, flexion);  % Moment arm
        ECRB.computeActuation(state);        % Force
        ECRB_force=ECRB.getTendonForce(state); % This is the total force which is the sum of active and passive and takes account of the pennation angle (if any)
        ECRB_active=ECRB.getActiveFiberForce(state); % Active force
        ECRB_passive=ECRB.getPassiveFiberForce(state); % Passive force
        
        % Brachioradialis
        % Moment arm:
        br1_ma=br1.computeMomentArm(state, flexion); % Br1
        br2_ma=br2.computeMomentArm(state, flexion); % Br2
        br3_ma=br3.computeMomentArm(state, flexion); % Br3
        
        % Force:
        br1_force=br1.getTendonForce(state);
        br1_active=br1.getActiveFiberForce(state);
        br1_passive=br1.getPassiveFiberForce(state);

        br2_force=br2.getTendonForce(state); % Br2
        br2_active=br2.getActiveFiberForce(state);
        br2_passive=br2.getPassiveFiberForce(state);

        br3_force=br3.getTendonForce(state); % Br3
        br3_active=br3.getActiveFiberForce(state);
        br3_passive=br3.getPassiveFiberForce(state);

        br_force=br1_force+br2_force+br3_force; % Total force sum
        br_active=br1_active+br2_active+br3_active; % Active force sum
        br_passive=br1_passive+br2_passive+br3_passive; % Passive force sum
        
        % Calculate moment 
        ECRB_moment=ECRB_ma*ECRB_force;

        % Active moment
        ecrb_moment_a=ECRB_ma*ECRB_active;
        
        % Passive moment
        ecrb_moment_p=ECRB_ma*ECRB_passive;

        % ALL EXTENSORS
        % MOMENT ARM
        ecu_ma=ECU.computeMomentArm(state, flexion); 
        edci_ma=EDCI.computeMomentArm(state, flexion); 
        edcl_ma=EDCL.computeMomentArm(state, flexion); 
        edcm_ma=EDCM.computeMomentArm(state, flexion); 
        edcr_ma=EDCR.computeMomentArm(state, flexion); 
        edm_ma=EDM.computeMomentArm(state, flexion); 
        eip_ma=EIP.computeMomentArm(state, flexion); 
        epl_ma=EPL.computeMomentArm(state, flexion); 
        % FORCE
        ecu_force=ECU.getTendonForce(state);
        edci_force=EDCI.getTendonForce(state);
        edcl_force=EDCL.getTendonForce(state);
        edcm_force=EDCM.getTendonForce(state);
        edcr_force=EDCR.getTendonForce(state);
        edm_force=EDM.getTendonForce(state);
        eip_force=EIP.getTendonForce(state);
        epl_force=EPL.getTendonForce(state);
        % MOMENT
        ecu_m=ecu_force*ecu_ma;
        edci_m=edci_force*edci_ma;
        edcl_m=edcl_ma*edcl_force;
        edcm_m=edcm_ma*edcm_force;
        edcr_m=edcr_ma*edcm_force;
        edm_m=edm_ma*edm_force;
        eip_m=eip_ma*eip_force;
        epl_m=epl_ma*epl_force;
        % TOTAL MOMENT
        total_moment=ecu_m+edci_m+edcl_m+edcm_m+edcr_m+edm_m+eip_m+epl_m+ECRB_moment;      
        
        % Store solutions
        position=f*180/pi; % Position in degrees
        ecrb_f=[ecrb_f;ECRB_ma,ECRB_active,ECRB_passive,ECRB_force,ECRB_moment,position,]; %ECRB results
        moment_f=[moment_f;ecrb_moment_a, ecrb_moment_p,ECRB_moment];
        total_m=[total_m;total_moment];

    end
    
% --------- ----------  RESET MODEL --------- ---------- --------- --------
    
    flexion.setValue(state, 0);

     for i=0:nMuscles-1
        muscles.get(i).setActivation(state,1.0)
        afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
     end

    osimModel.equilibrateMuscles(state); 
    
% --------- ----------  ---------- --------- ---------- --------- ---------


  % ELBOW FLEXION: step through 100 points in ROM:
    for e=linspace(rMinE,rMaxE,100)
        elbow.setValue(state, e);

        for i=0:nMuscles-1 
             muscles.get(i).setActivation(state,1.0)
             afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
        end

        osimModel.equilibrateMuscles(state); 

        %-------------------------------------------------------------ECRB:
        ECRB_ma_w=ECRB.computeMomentArm(state, flexion); % Wrist flexion moment arm
        ECRB_ma_e=ECRB.computeMomentArm(state, elbow); %Elbow flexion moment arm
        ECRB.computeActuation(state); % Force
        ECRB_force=ECRB.getTendonForce(state); % This is the total force which is the sum of active and passive and takes account of the pennation angle (if any)
        ECRB_active=ECRB.getActiveFiberForce(state); % Active force
        ECRB_passive=ECRB.getPassiveFiberForce(state); % Passive force

        %--------------------------------------------------Brachioradialis:

        % Moment arm:
        % Wrist Flexion:
        br1_ma_w=br1.computeMomentArm(state, flexion); % Br1 
        br2_ma_w=br2.computeMomentArm(state, flexion); % Br2
        br3_ma_w=br3.computeMomentArm(state, flexion); % Br3
        % Elbow Flexion
        br1_ma_e=br1.computeMomentArm(state, elbow); % Br1
        br2_ma_e=br2.computeMomentArm(state, elbow); % Br2
        br3_ma_e=br3.computeMomentArm(state, elbow); % Br3

        % Force:
        br1_force=br1.getTendonForce(state);   %---------Br1
        br1_active=br1.getActiveFiberForce(state); 
        br1_passive=br1.getPassiveFiberForce(state);
        br2_force=br2.getTendonForce(state);   %---------Br2
        br2_active=br2.getActiveFiberForce(state);
        br2_passive=br2.getPassiveFiberForce(state);
        br3_force=br3.getTendonForce(state);   %---------Br3
        br3_active=br3.getActiveFiberForce(state);
        br3_passive=br3.getPassiveFiberForce(state);
        br_force=br1_force+br2_force+br3_force; % Total force sum
        br_active=br1_active+br2_active+br3_active; % Active force sum
        br_passive=br1_passive+br2_passive+br3_passive; % Passive force sum

        % Calculate moment 
        ECRB_moment_w=ECRB_ma_w*ECRB_force; 
        br1_moment_w=br1_ma_w*br1_force; 
        br2_moment_w=br2_ma_w*br2_force; 
        br3_moment_w=br3_ma_w*br3_force; 
        br_moment_w=br1_moment_w+br2_moment_w+br3_moment_w; % Moment sum
        br1_moment_e=br1_ma_e*br1_force;
        br2_moment_e=br2_ma_e*br2_force;
        br3_moment_e=br3_ma_e*br3_force;
        br_moment_e=br1_moment_e+br2_moment_e+br3_moment_e;
        
         % Active moment
        br1_moment_a=br1_ma_e*br1_active;
        br2_moment_a=br2_ma_e*br2_active;
        br3_moment_a=br3_ma_e*br3_active;
        br_moment_a=br1_moment_a+br2_moment_a+br3_moment_a;
        
        % Passive moment
        br1_moment_p=br1_ma_e*br1_passive;
        br2_moment_p=br2_ma_e*br2_passive;
        br3_moment_p=br3_ma_e*br3_passive;
        br_moment_p=br1_moment_p+br2_moment_p+br3_moment_p;
       
        % ALL EXTENSORS
        
        % MOMENT ARM
        brach_1_ma=brachialis_1.computeMomentArm(state, elbow); 
        brach_2_ma=brachialis_2.computeMomentArm(state, elbow); 
        brach_3_ma=brachialis_3.computeMomentArm(state, elbow); 
        brach_4_ma=brachialis_4.computeMomentArm(state, elbow); 
        brach_5_ma=brachialis_5.computeMomentArm(state, elbow); 
        brach_6_ma=brachialis_6.computeMomentArm(state, elbow); 
        brach_7_ma=brachialis_7.computeMomentArm(state, elbow); 
        bic_1_ma=bic_b_1.computeMomentArm(state, elbow);         
        bic_2_ma=bic_b_2.computeMomentArm(state, elbow); 
        % FORCE
        brach_f_1=brachialis_1.getTendonForce(state);
        brach_f_2=brachialis_2.getTendonForce(state);
        brach_f_3=brachialis_3.getTendonForce(state);
        brach_f_4=brachialis_4.getTendonForce(state);
        brach_f_5=brachialis_5.getTendonForce(state);
        brach_f_6=brachialis_6.getTendonForce(state);
        brach_f_7=brachialis_7.getTendonForce(state);
        bic_f_1=bic_b_1.getTendonForce(state);
        bic_f_2=bic_b_2.getTendonForce(state);
        % MOMENT
        brach_m_1=brach_1_ma*brach_f_1;
        brach_m_2=brach_2_ma*brach_f_2;
        brach_m_3=brach_3_ma*brach_f_3;
        brach_m_4=brach_4_ma*brach_f_4;
        brach_m_5=brach_5_ma*brach_f_5;
        brach_m_6=brach_6_ma*brach_f_6;
        brach_m_7=brach_7_ma*brach_f_7;
        bic_m_1=bic_f_1*bic_1_ma;
        bic_m_2=bic_f_2*bic_2_ma;
        
        el_flex=[el_flex;(brach_m_1+brach_m_2+brach_m_3+brach_m_4+brach_m_5+brach_m_6+brach_m_7+bic_m_1+bic_m_2+br_moment_e)];

        % Store solutions
        position=e*180/pi; % Position in degrees
        ECRB_e=[ECRB_e;ECRB_ma_w,ECRB_active,ECRB_passive,ECRB_force,ECRB_moment_w,position,ECRB_ma_e]; %ECRB results
        br_e=[br_e;br1_ma_w,br2_ma_w,br3_ma_w,br1_force,br2_force,br3_force,br1_moment_w,br2_moment_w,br3_moment_w,br_active,br_passive,br_force,br_moment_w,position,br_moment_e,br1_ma_e]; %Br results
        moment_e=[moment_e;br_moment_a,br_moment_p,br_moment_e];
    end
    
% -------------- ----------  RESET MODEL --------- ---------- ----------  
        flexion.setValue(state, 0);

         for i=0:nMuscles-1
            muscles.get(i).setActivation(state,1.0)
            afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
         end

        osimModel.equilibrateMuscles(state); 

        elbow.setValue(state, 30);

         for i=0:nMuscles-1
            muscles.get(i).setActivation(state,1.0)
            afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
         end

% -------------- ----------  ---------- --------- ---------- ---------

        osimModel.equilibrateMuscles(state); 
    
         % DEVIATION: step through 100 points in ROM: 
     for d=linspace(rMinDev,rMaxDev,100)

         deviation.setValue(state, d);         

         for i=0:nMuscles-1 
             muscles.get(i).setActivation(state,1.0)
             afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
         end

        osimModel.equilibrateMuscles(state); 

        %-----------------------------------------------------ECRB:
        ECRB_ma=ECRB.computeMomentArm(state, deviation);  % Moment arm
        ECRB.computeActuation(state);        % Force
        ECRB_force=ECRB.getTendonForce(state); % This is the total force which is the sum of active and passive and takes account of the pennation angle (if any)
        ECRB_active=ECRB.getActiveFiberForce(state); % Active force
        ECRB_passive=ECRB.getPassiveFiberForce(state); % Passive force

        %------------------------------------------Brachioradialis:
        % Moment arm:
        br1_ma=br1.computeMomentArm(state, deviation); % Br1
        br2_ma=br2.computeMomentArm(state, deviation); % Br2
        br3_ma=br3.computeMomentArm(state, deviation); % Br3

        % Force:
        br1_force=br1.getTendonForce(state);
        br1_active=br1.getActiveFiberForce(state);
        br1_passive=br1.getPassiveFiberForce(state);

        br2_force=br2.getTendonForce(state); % Br2
        br2_active=br2.getActiveFiberForce(state);
        br2_passive=br2.getPassiveFiberForce(state);

        br3_force=br3.getTendonForce(state); % Br3
        br3_active=br3.getActiveFiberForce(state);
        br3_passive=br3.getPassiveFiberForce(state);

        br_force=br1_force+br2_force+br3_force; % Total force sum
        br_active=br1_active+br2_active+br3_active; % Active force sum
        br_passive=br1_passive+br2_passive+br3_passive; % Passive force sum
        
        % Calculate moment 
        ECRB_moment=ECRB_ma*ECRB_force;

        % Active moment
        ecrb_moment_a=ECRB_ma*ECRB_active;
        
        % Passive moment
        ecrb_moment_p=ECRB_ma*ECRB_passive;
        
        % Calculate moment 
        ECRB_moment=ECRB_ma*ECRB_force;
        
        
        
        % Store solutions
        position=d*180/pi; % Position in degrees
        ecrb_d=[ecrb_d;ECRB_ma,ECRB_active,ECRB_passive,ECRB_force,ECRB_moment,position,]; %ECRB results
        moment_d=[moment_d;ecrb_moment_a,ecrb_moment_p,ECRB_moment];
        end

     
%----------ECRB values----------------------------------------------------

% Force
ecrb_pre_f_active=ecrb_f(:,2);
ecrb_pre_f_passive=ecrb_f(:,3);
ecrb_pre_f=ecrb_f(:,4);
% Moment
ecrb_pre_m=ecrb_f(:,5);
ecrb_pre_m_a=moment_f(:,1);
ecrb_pre_m_p=moment_f(:,2);
save('ecrb_pre_m.mat','ecrb_pre_m');
% Deviation
dev_pre_m=ecrb_d(:,5);
dev_pre_m_a=moment_f(:,1);
dev_pre_m_p=moment_f(:,2);


%-------BR Values---------------------------------------------------------
% Force
br_pre_f_active=br_e(:,10);
br_pre_f_passive=br_e(:,11);
br_pre_f_total=br_e(:,12);

% Moment
el_pre_m_a=moment_e(:,1);
el_pre_m_p=moment_e(:,2);
el_pre_m=moment_e(:,3);

%--------------------------------------------------------------------------

wrist_deviation=ecrb_d(1:100,6);        
wrist_flexion=ecrb_f(1:100,6);
wrist_extension=wrist_flexion*-1;
elbow_flexion=br_e(1:100,14);

%-----------WRIST FLEXION Investigating change with wrist flexion----------

%---       ---      Figure 1 (force and moment)     ---       ---       ---
% 
% figure('name','Force/ Moment of Brachioradialis and ECRB')
% tiledlayout('flow')
% 
% % ECRB force with active/passive components
% nexttile
% plot(wrist_extension,ecrb_f(1:100,2)); hold on
% plot(wrist_extension,ecrb_f(1:100,3)); hold on
% plot(wrist_extension,ecrb_f(1:100,4));
% title('ECRB Force Pre-transfer')
% xlabel('Wrist Extension (°)')
% ylabel('Force (N)')
% legend('Active','Passive','Total')
% 
% 
% % Br force with active/passive components
% nexttile
% plot(elbow_flexion,br_e(1:100,10)); hold on
% plot(elbow_flexion,br_e(1:100,11)); hold on
% plot(elbow_flexion,br_e(1:100,12));
% title('Brachioradialis Force Pre-transfer')
% xlabel('Elbow Flexion (°)')
% ylabel('Force (N)')
% legend('Active','Passive','Total')
% 
% ECRB Wrist Flexion Moment 
nexttile
plot(wrist_extension,(-1*ecrb_f(1:100,5))); 
title('ECRB Wrist Flexion Moment Pre-transfer')
xlabel('Wrist Extension (°)')
ylabel('Wrist Extension Moment (Nm)')

% 
% % BR Elbow flexion moment pre
% nexttile
% plot(elbow_flexion,br_e(1:100,15)); hold on
% plot(elbow_flexion,el_pre_m_a); hold on
% plot(elbow_flexion,el_pre_m_p); hold on
% title('Brachioradialis Elbow Flexion Moment Pre-transfer')
% legend('total','active','passive')
% xlabel('Elbow Flexion (°)')
% ylabel('Elbow Flexion Moment (Nm)')
% 
% % ECRB Wrist deviation moment pre
% nexttile
% plot(wrist_deviation,ecrb_d(1:100,5))
% title('ECRB Deviation Moment Pre-tranfer')
% xlabel('Wrist Deviation (°)')
% ylabel('Deviation Moment (Nm)')
% 
nexttile
plot(wrist_flexion,total_m)

nexttile
plot(elbow_flexion,el_flex); hold on
plot(elbow_flexion,el_pre_m);
