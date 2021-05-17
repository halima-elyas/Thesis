clc; clear; close all;


%--------------------------------------------------------------------------
set(0,'DefaultaxesFontSize', 12) 
set(0,'DefaultTextFontSize', 12) 
set(0,'DefaultaxesFontName', 'Times new Roman') 
set(0,'DefaultlegendFontName', 'Times new Roman')
set(0,'defaultAxesXGrid','on') 
set(0,'defaultAxesYGrid','on') 
%--------------------------------------------------------------------------

%PRE TRANSFER VALUES
total_m_pre=-1*[-9.24863853809704;-9.75693730422713;-9.84448781989014;-9.92580083681134;-10.0020015838505;-9.65571814830926;-9.73360637212079;-10.1968187978954;-9.89596421260673;-10.2978798193283;-10.4861264515341;-10.7956112693509;-10.6896260175977;-11.1341130729745;-11.1875283525669;-11.2310212538462;-11.2658014619463;-10.7784069832951;-10.7691250339536;-11.3194300748019;-11.3170471209857;-11.3066893475203;-11.2857430929011;-11.2550009279811;-11.2149426915510;-11.1656894886571;-11.1104950260861;-11.0459155066040;-10.9737651986360;-10.8936494865564;-10.8096482322604;-10.7184044317659;-10.6225568513267;-10.5215212991763;-10.4159073968507;-10.3065639487952;-10.1940316027002;-10.0771171019601;-9.95723666696008;-9.83602869768004;-9.71389065044409;-9.59056854115371;-9.46668156386375;-9.34187217083773;-9.21697891505582;-9.09148809662550;-8.96526493288022;-8.83920354613280;-8.71364049588307;-8.58800875836645;-8.46296844644375;-8.33830751278786;-8.21358969409352;-8.08831481663765;-7.96594897094912;-7.84320344271156;-7.72171042890420;-7.60087460517829;-7.48035021213822;-7.36031951861195;-7.24135246281172;-7.12370414757682;-7.00696633456968;-6.89033458508277;-6.77526140273424;-6.66078699047535;-6.54666886065847;-6.43327794059305;-6.32094512064029;-6.20948065313676;-6.09800287150780;-5.98760199691633;-5.87794951713125;-5.76870619913750;-5.66032334693255;-5.55237255808185;-5.44470100214504;-5.33819892202583;-5.23142988243478;-5.12650546320448;-5.02159045800678;-4.91763709322671;-4.81412159224699;-4.71205272778260;-4.61016161775287;-4.50885904636476;-4.40883022591137;-4.30894101542808;-4.21017230748343;-4.11171386349250;-4.01344803828957;-3.91610456511752;-3.81996457821819;-3.72392030396402;-3.62852318209892;-3.53356026676099;-3.43883114577624;-3.34544621289569;-3.25238076946121;-3.15965146353290];
ecrb_m_pre=-1*[-2.81004587417053;-2.79589246116414;-2.78079559763814;-2.76475450099692;-2.74777125993052;-2.72987132675708;-2.71108468671689;-2.69144047743479;-2.67096689249862;-2.64957605303926;-2.62781780189408;-2.60498730431922;-2.58143290718065;-2.55718028640176;-2.53225484526034;-2.50668311325779;-2.48048963233973;-2.45387896875864;-2.42650412311492;-2.39857877084111;-2.37016876050841;-2.34121886342412;-2.31179481735985;-2.28192813208635;-2.25162206893159;-2.22101903748404;-2.19056346850747;-2.16044787390367;-2.13030876569863;-2.10034073081870;-2.07064812650241;-2.04100719800840;-2.01155698491819;-1.98230206552462;-1.95324663233035;-1.92439448093070;-1.89589084017013;-1.86735292491738;-1.83923143702910;-1.81120898700041;-1.78340362310073;-1.75585019773119;-1.72849884027023;-1.70138762114943;-1.67451846974236;-1.64789458300804;-1.62153132265749;-1.59543490802403;-1.56961236469730;-1.54410598725404;-1.51885238361227;-1.49395907807268;-1.46929397590466;-1.44492100159890;-1.42086936010483;-1.39713038385253;-1.37370763967886;-1.35066160373357;-1.32787624051211;-1.30541537556931;-1.28316410156805;-1.26130740262656;-1.23977254939212;-1.21855995768441;-1.19766969785089;-1.17710150847410;-1.15685481011935;-1.13696305131359;-1.11735324840450;-1.09804936528478;-1.07907508516166;-1.06041546387594;-1.04206852840312;-1.02403561360473;-1.00632616221884;-0.988949541538658;-0.971913618034659;-0.955224810960339;-0.938888150425718;-0.922907339182503;-0.907284817383785;-0.892021829610119;-0.877118493491406;-0.862438474191067;-0.848206579017085;-0.834421890209055;-0.820910470080975;-0.807741238295830;-0.794908795654501;-0.782407111404160;-0.770229589299407;-0.758369131331155;-0.746818198820803;-0.735568870053198;-0.724612896667844;-0.713941758333231;-0.703546711397899;-0.693369224812036;-0.683512516275329;-0.673928362650886];
br_m_pre=[-1.00599942380051;-0.855629600049860;-0.707233825300961;-0.560499815913360;-0.415108456716555;-0.270740590100443;-0.127078428548780;0.0161941775960124;0.159391559622745;0.302826110578037;0.446808197670290;0.591646078687227;0.737645793791630;0.885110928321029;1.03433633887786;1.18561116391068;1.33921891075819;1.49543533296048;1.65452719973941;1.81675075197337;1.98234999429622;2.15155475872654;2.32457850476102;2.50157633209974;2.68284001663743;2.86839879495208;3.05841290709951;3.25312946632143;3.45224586147346;3.65593963947515;3.86416302665745;4.07685850616784;4.29386008955884;4.51476105723609;4.73903292799202;4.96631990083375;5.19620803056848;5.42822039643221;5.66181461777057;5.89673362064952;6.13151568950062;6.36641293697996;6.59934372458871;6.83020128636720;7.05819831806937;7.28253528715801;7.50238836942769;7.71692057581780;7.92536139538873;8.12699369786108;8.32087490027650;8.50670989198518;8.68451942106332;8.85311467913293;9.01213081607525;9.16187748502785;9.30199508931938;9.43188055506630;9.55167160980587;9.66140741105402;9.76121408261149;9.85144115687575;9.93259672123712;10.0049488317308;10.0702214632880;10.1282201949493;10.1802242566503;10.2265667677539;10.2669442909755;10.3002250559721;10.3251049322989;10.3402190089549;10.3441551585494;10.3354730877578;10.3127250781874;10.2744719846864;10.2189049743987;10.1435446945165;10.0455581079950;9.92193094476860;9.76977912820050;9.58636522329920;9.36876763877171;9.11412722139365;8.82013888483266;8.48784488992992;8.12104822660526;7.72361769699645;7.29922622304820;6.85160441860466;6.38601931081558;5.91365411405242;5.44503441801396;4.98847386469645;4.54571885931818;4.11945932924778;3.71238904330591;3.33124793347746;2.97800432101285;2.65831352600200];
dev_m_pre=-1*[-3.47480581036268;-3.42761361557070;-3.38079687536398;-3.33556677768435;-3.29110949796436;-3.24772276369384;-3.20491835820643;-3.16314216261262;-3.12232276696010;-3.08246277085916;-3.04356137843216;-3.00561466541311;-2.96861583332622;-2.93255545296493;-2.89742170764333;-2.86369987967139;-2.83037160301644;-2.79790751099107;-2.76633197929726;-2.73558064252533;-2.70564613613297;-2.67650560321965;-2.64813534904743;-2.62051212489886;-2.59339531771829;-2.56719666343598;-2.54167766100824;-2.51681330039908;-2.49257451367133;-2.46893003550263;-2.44540579298895;-2.42286563933599;-2.40082917502330;-2.37926729232936;-2.35815150576089;-2.33726207155842;-2.31692322267878;-2.29694781235051;-2.27731045292992;-2.25798647244935;-2.23895191741196;-2.22018355372226;-2.20165886592263;-2.18335605490072;-2.16525403421464;-2.14708194701645;-2.12921112789874;-2.11160416788179;-2.09396379550323;-2.07659612019224;-2.05931725958384;-2.04195075790921;-2.02481053650602;-2.00771205558152;-1.99064092783490;-1.97358338052169;-1.95652624347716;-1.93945693691869;-1.92236345907919;-1.90523437370657;-1.88805879747284;-1.87082638732658;-1.85340550947274;-1.83604475707030;-1.81859892431356;-1.80105970083028;-1.78341925369952;-1.76557128201460;-1.74772989730937;-1.72976587100322;-1.71167316013742;-1.69344613855059;-1.67507958513605;-1.65656867226600;-1.63790895439645;-1.61909635685418;-1.60012716481847;-1.58099029420783;-1.56166566558500;-1.54217751227574;-1.52252343383543;-1.50260994282470;-1.48261012810258;-1.46243925151194;-1.44209604210740;-1.42157949157125;-1.40088884642788;-1.38051183390883;-1.35950952719774;-1.33833303677546;-1.31692871897862;-1.29495413454763;-1.27322799332357;-1.25068218146329;-1.22857801060145;-1.20630291747604;-1.18385835455020;-1.16084079181620;-1.13804446058129;-1.11508476604793];
el_m_pre=[10.9332624515142;11.4173306027104;11.9115256386288;12.4169638573808;12.9353405547877;13.4674993166595;14.1254730397239;14.8546956103151;15.6043740285386;16.3783974184623;17.3756226125622;18.4453250163210;19.5552833696724;20.6971716888081;21.8716236882223;23.0921344224025;24.3613107124802;25.6653441074190;27.0025496235314;28.3710147764093;29.7700253619971;31.1969220466130;32.6488732584178;34.1541760794933;35.6832692551287;37.2278333040895;38.7822253094062;40.3418354908069;41.8981720786415;43.4448763416499;44.9764076615947;46.4873073018267;47.9725488398696;49.4297853487393;50.8527408173514;52.2380524985116;53.5816463186212;54.8793956716475;56.1288624019562;57.3265790792758;58.4682393007906;59.5527362756332;60.5757419372841;61.5357403057269;62.4302397383204;63.2575752656317;64.0147717008912;64.6994711031401;65.3127515947655;65.8538235936878;66.3222386033240;66.7152189247557;67.0320193530383;67.2700300766811;67.4272930970644;67.5036514082701;67.4966894304813;67.4046455351550;67.2267042643070;66.9623341729768;66.6115512089947;66.1747609320037;65.6536769714994;65.0522721305506;64.3714710211389;63.6091122015412;62.7646569493686;61.8374678490622;60.8252062972127;59.7248910584594;58.5334778704357;57.2412343525542;55.8503071672070;54.3636260938611;52.7839025918053;51.1168147747335;49.3697085805302;47.5509410913706;45.6707581414194;43.7396329117575;41.7722330643401;39.7833979033643;37.7767897596576;35.7470392156268;33.6916244492270;31.6183860766864;29.5425821517910;27.4833478675041;25.4606522718628;23.4945628434034;21.6063761262259;19.8201876548024;18.1497257462323;16.5927353373137;15.1400676898272;13.7842784882598;12.5169334144101;11.3386118730904;10.2479191706721;9.24709032843888];

%--------------------------------------------------------------------------

import org.opensim.modeling.*

% Read in osim model
osimModel=Model('arm_post_transfer_29_30.osim');    % create model object by loading the model from file
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
Br1geopath=br1.getGeometryPath();
 
% Get a handle to the ECRB
ECRB = osimModel.getMuscles().get("ECRB");
ECRBgeopath=ECRB.getGeometryPath();
ECRBpathpoint=ECRBgeopath.getPathPointSet().get(0);

% Get a handle to rest of wrist moving muscles
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

%--------------------------------------------------------------------------

% Backup tendon slack length 
backup_brach1_TendonSlackLength=br1.getTendonSlackLength();     
backup_brach2_TendonSlackLength=br2.getTendonSlackLength();     
backup_brach3_TendonSlackLength=br3.getTendonSlackLength();     

% Assign slack length variables
br1_sl=backup_brach1_TendonSlackLength;
br2_sl=backup_brach2_TendonSlackLength;
br3_sl=backup_brach3_TendonSlackLength;

%-----------------------------Solution------------------------------------%

% Define solution matrices 
ecrb_f=[];
br_slack_f=[];
br_f_ap=[];

ecrb_d=[];
br_slack_d=[];
br_d=[];
br_d_ap=[];
br_ma_f=[];
ECRB_el=[];
br_slack_el=[];
br_el=[];
br_el_ap=[];
br_moment_f=[];
br_moment_el=[];
br_moment_d=[];
total_moment=[];
el_flex=[];

for i=0:nMuscles-1
muscles.get(i).setActivation(state,1.0)
afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i)); 
end

for s=1:6
    
    br1.setTendonSlackLength(br1_sl);
    br1_sl=br1.getTendonSlackLength();
    br2.setTendonSlackLength(br2_sl);
    br2_sl=br2.getTendonSlackLength();
    br3.setTendonSlackLength(br3_sl);
    br3_sl=br3.getTendonSlackLength();

    osimModel.initSystem(); % Re-initialise system when a parameter is changed

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
        ECRB_ma=ECRB.computeMomentArm(state, flexion)*-1;  % Moment arm
        ECRB.computeActuation(state);        % Force
        ECRB_force=ECRB.getTendonForce(state); % This is the total force which is the sum of active and passive and takes account of the pennation angle (if any)
        ECRB_active=ECRB.getActiveFiberForce(state); % Active force
        ECRB_passive=ECRB.getPassiveFiberForce(state); % Passive force

        % Brachioradialis
        % Moment arm:
        br1_ma=br1.computeMomentArm(state, flexion)*-1; % Br1
        br2_ma=br2.computeMomentArm(state, flexion)*-1; % Br2
        br3_ma=br3.computeMomentArm(state, flexion)*-1; % Br3

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
        br1_moment=br1_ma*br1_force;
        br2_moment=br2_ma*br2_force;
        br3_moment=br3_ma*br3_force;
        br_moment=br1_moment+br2_moment+br3_moment; % Moment sum

          %Active moment
        br1_moment_a=br1_ma*br1_active;
        br2_moment_a=br2_ma*br2_active;
        br3_moment_a=br3_ma*br3_active;
        br_moment_a=br1_moment_a+br2_moment_a+br3_moment_a;
        
        %Passive moment
        br1_moment_p=br1_ma*br1_passive;
        br2_moment_p=br2_ma*br2_passive;
        br3_moment_p=br3_ma*br3_passive;
        br_moment_p=br1_moment_p+br2_moment_p+br3_moment_p;
        
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
        total_moment=[total_moment;ecu_m+edci_m+edcl_m+edcm_m+edcr_m+edm_m+eip_m+epl_m+br_moment];
    
        if br_moment_a>=br_moment_p
            ROM=1;
        else
            ROM=0;
        end
        
        % Get fiber length
        
        % Store solutions
        position=f*180/pi; % Position in degrees
        ecrb_f=[ecrb_f;ECRB_ma,ECRB_active,ECRB_passive,ECRB_force,ECRB_moment,position,]; %ECRB results
        br_slack_f=[br_slack_f;br1_ma,br2_ma,br3_ma,br1_force,br2_force,br3_force,br1_moment,br2_moment,br3_moment,br_active,br_passive,br_force,br_moment,position,br1_sl,br2_sl,br3_sl]; %All br results
        br_f=[br_f;br1_ma,br2_ma,br3_ma,br_force,br_moment,position]; % Relevant BR values (moment arm, total force and total moment)
        br_f_ap=[br_f_ap;br1_active,br1_passive,br1_force,br2_active,br2_passive,br2_force,br3_active,br3_passive,br3_force,br_active,br_passive,br_force,position]; %BR passive/ active split
        br_moment_f=[br_moment_f;br_moment_a,br_moment_p,br_moment,ROM];

    end

flexion.setValue(state, 0);

     for i=0:nMuscles-1
        muscles.get(i).setActivation(state,1.0)
        afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
     end

    osimModel.equilibrateMuscles(state); 

     %ELBOW FLEXION: step through 100 points in ROM:
    for e=linspace(rMinE,rMaxE,100)
        elbow.setValue(state, e);

        for i=0:nMuscles-1 
             muscles.get(i).setActivation(state,1.0)
             afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
        end

        osimModel.equilibrateMuscles(state); 

        %-------------------------------------------------------------ECRB:
        ECRB_ma_w=ECRB.computeMomentArm(state, flexion)*-1; % Wrist Extension moment arm
        ECRB_ma_e=ECRB.computeMomentArm(state, elbow); %Elbow flexion moment arm
        ECRB.computeActuation(state); % Force
        ECRB_force=ECRB.getTendonForce(state); % This is the total force which is the sum of active and passive and takes account of the pennation angle (if any)
        ECRB_active=ECRB.getActiveFiberForce(state); % Active force
        ECRB_passive=ECRB.getPassiveFiberForce(state); % Passive force

        %--------------------------------------------------Brachioradialis:

        % Moment arm:
        % Wrist Extension:
        br1_ma_w=br1.computeMomentArm(state, flexion)*-1; % Br1 
        br2_ma_w=br2.computeMomentArm(state, flexion)*-1; % Br2
        br3_ma_w=br3.computeMomentArm(state, flexion)*-1; % Br3
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

        %Active moment
        br1_moment_a=br1_ma_w*br1_active;
        br2_moment_a=br2_ma_w*br2_active;
        br3_moment_a=br3_ma_w*br3_active;
        br_moment_a=br1_moment_a+br2_moment_a+br3_moment_a;
        
        %Passive moment
        br1_moment_p=br1_ma_w*br1_passive;
        br2_moment_p=br2_ma_w*br2_passive;
        br3_moment_p=br3_ma_w*br3_passive;
        br_moment_p=br1_moment_p+br2_moment_p+br3_moment_p;
       
        if br_moment_a>=br_moment_p
            ROM=1;
        else
            ROM=0;
        end

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

        el_flex=[el_flex;brach_m_1+brach_m_2+brach_m_3+brach_m_4+brach_m_5+brach_m_6+brach_m_7+bic_m_1+bic_m_2+br_moment_e];
            
        % Store solutions
        position=e*180/pi; % Position in degrees
        ECRB_el=[ECRB_el;ECRB_ma_w,ECRB_active,ECRB_passive,ECRB_force,ECRB_moment_w,position,ECRB_ma_e]; %ECRB results
        br_slack_el=[br_slack_el;br1_ma_w,br2_ma_w,br3_ma_w,br1_force,br2_force,br3_force,br1_moment_w,br2_moment_w,br3_moment_w,br_active,br_passive,br_force,br_moment_w,position,br_moment_e,br1_ma_e]; %Br results
        br_el=[br_el;br_force,br_moment_w,br_moment_e];
        br_el_ap=[br_el_ap;br_active,br_passive,br_force,position];
        br_moment_el=[br_moment_el;br_moment_a,br_moment_p,br_moment,ROM]; % Active/ passive comp. of wrist flex moment

    end
    
    flexion.setValue(state, 0);

     for i=0:nMuscles-1
        muscles.get(i).setActivation(state,1.0)
        afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
     end

    osimModel.equilibrateMuscles(state); 

    elbow.setValue(state, 0.523599);

     for i=0:nMuscles-1
        muscles.get(i).setActivation(state,1.0)
        afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
     end

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
        ECRB_ma=ECRB.computeMomentArm(state, deviation)*-1;  % Moment arm
        ECRB.computeActuation(state);        % Force
        ECRB_force=ECRB.getTendonForce(state); % This is the total force which is the sum of active and passive and takes account of the pennation angle (if any)
        ECRB_active=ECRB.getActiveFiberForce(state); % Active force
        ECRB_passive=ECRB.getPassiveFiberForce(state); % Passive force

        %------------------------------------------Brachioradialis:
        % Moment arm:
        br1_ma=br1.computeMomentArm(state, deviation)*-1; % Br1
        br2_ma=br2.computeMomentArm(state, deviation)*-1; % Br2
        br3_ma=br3.computeMomentArm(state, deviation)*-1; % Br3

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
        br1_moment=br1_ma*br1_force;
        br2_moment=br2_ma*br2_force;
        br3_moment=br3_ma*br3_force;
        br_moment=br1_moment+br2_moment+br3_moment; % Moment sum

        % Store solutions
        position=d*180/pi; % Position in degrees
        ecrb_d=[ecrb_d;ECRB_ma,ECRB_active,ECRB_passive,ECRB_force,ECRB_moment,position,]; %ECRB results
        br_slack_d=[br_slack_d;br1_ma,br2_ma,br3_ma,br1_force,br2_force,br3_force,br1_moment,br2_moment,br3_moment,br_active,br_passive,br_force,br_moment,position]; %All br results
        br_d=[br_d;br1_ma,br2_ma,br3_ma,br_force,br_moment,position]; % Relevant BR values (moment arm, total force and total moment)
        br_d_ap=[br_d_ap;br1_active,br1_passive,br1_force,br2_active,br2_passive,br2_force,br3_active,br3_passive,br3_force,br_active,br_passive,br_force,position]; %BR passive/ active split
        br_moment_d=[br_moment_d;br_moment_a,br_moment_p,br_moment,ROM];

     end
     
        % Update slack length 
        br1_sl=br1_sl+0.02
        br2_sl=br2_sl+0.02;
        br3_sl=br3_sl+0.02; 
     
end

wrist_deviation=-1*br_d(1:100,6);
wrist_flexion=-1*br_f(1:100,6);
wrist_extension=wrist_flexion;
elbow_flexion=br_el_ap(1:100,4);

%--------------------------------------------------------------------------
%---------------------------Wrist Extension----------------------------------
%--------------------------------------------------------------------------

%---       ---      Figure 1 (force and moment)     ---       ---       ---

figure('name','Force/ Moment of Brachioradialis and ECRB')
tiledlayout('flow')

% ECRB force at 30 degrees elbow flex
nexttile
plot(wrist_extension,ecrb_f(1:100,4)); hold on
plot(wrist_extension,ecrb_f(101:200,4)); hold on
plot(wrist_extension,ecrb_f(201:300,4)); hold on
plot(wrist_extension,ecrb_f(301:400,4)); hold on
plot(wrist_extension,ecrb_f(401:500,4));
title('ECRB at 30 degrees Wrist Extension')
xlabel('Wrist Extension (°)')
ylabel('Force (N)') 
legend('Length 5','Length 6','Length 7','Length 8','Length 9')

% ECRB Wrist Extension Moment at 30 degrees elbow flex
nexttile
plot(wrist_extension,ecrb_f(1:100,5)); hold on
plot(wrist_extension,ecrb_f(101:200,5)); hold on
plot(wrist_extension,ecrb_f(201:300,5)); hold on
plot(wrist_extension,ecrb_f(301:400,5)); hold on
plot(wrist_extension,ecrb_f(401:500,5)); 
title('ECRB at 30 degrees Wrist Extension')
xlabel('Wrist Extension (°)')
ylabel('Wrist Extension Moment (Nm)')
legend('Length 5','Length 6','Length 7','Length 8','Length 9')

% ECRB Deviation Moment at 30 degrees elbow flex
nexttile
plot(wrist_deviation,ecrb_d(1:100,5)); hold on
plot(wrist_deviation,ecrb_d(101:200,5)); hold on
plot(wrist_deviation,ecrb_d(201:300,5)); hold on
plot(wrist_deviation,ecrb_d(301:400,5)); hold on
plot(wrist_deviation,ecrb_d(401:500,5)); 
title('ECRB at 30 degrees Wrist Extension')
xlabel('Wrist Deviation (°)')
ylabel('Wrist Deviation Moment (Nm)')
% legend('Length 5','Length 6','Length 7','Length 8','Length 9')

%BR Force at 30 deg
nexttile
plot(wrist_extension,br_f(1:100,4)); hold on
plot(wrist_extension,br_f(101:200,4)); hold on
plot(wrist_extension,br_f(201:300,4)); hold on
plot(wrist_extension,br_f(301:400,4)); hold on
plot(wrist_extension,br_f(401:500,4)); hold on
plot(wrist_extension,br_f(501:600,4)); 
xlabel('Wrist Extension (°)')
ylabel('Total Force (N)')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')

%BR Wrist Extension Moment at 30 degrees elbow flex
nexttile
plot(wrist_extension,br_f(1:100,5)); hold on
plot(wrist_extension,br_f(101:200,5)); hold on
plot(wrist_extension,br_f(201:300,5)); hold on
plot(wrist_extension,br_f(301:400,5)); hold on
plot(wrist_extension,br_f(401:500,5)); hold on
plot(wrist_extension,br_f(501:600,5)); hold on
plot(wrist_extension,ecrb_m_pre)
xlabel('Wrist Extension (°)')
ylabel('Total Wrist Extension Moment (Nm)')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6','ECRB unparalysed')

%Normalised Extension Moment
i1_normalised_w=br_f(1:100,5)./ecrb_m_pre;
i2_normalised_w=br_f(101:200,5)./ecrb_m_pre;
i3_normalised_w=br_f(201:300,5)./ecrb_m_pre;
i4_normalised_w=br_f(301:400,5)./ecrb_m_pre;
i5_normalised_w=br_f(401:500,5)./ecrb_m_pre;
i6_normalised_w=br_f(501:600,5)./ecrb_m_pre;

nexttile
plot(wrist_extension,(i1_normalised_w)); hold on
plot(wrist_extension,(i2_normalised_w)); hold on
plot(wrist_extension,(i3_normalised_w)); hold on
plot(wrist_extension,(i4_normalised_w)); hold on
plot(wrist_extension,(i5_normalised_w)); hold on
plot(wrist_extension,(i6_normalised_w)); hold on
plot(wrist_extension,ecrb_m_pre./ecrb_m_pre);
xlabel('Wrist Extension (°)');
ylabel('Normalised Wrist Extension Moment')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6','ECRB unparalysed')

%Active extension moment
nexttile
plot(wrist_extension,br_moment_f(1:100,1)); hold on
plot(wrist_extension,br_moment_f(101:200,1)); hold on
plot(wrist_extension,br_moment_f(201:300,1)); hold on
plot(wrist_extension,br_moment_f(301:400,1)); hold on
plot(wrist_extension,br_moment_f(401:500,1)); hold on
plot(wrist_extension,br_moment_f(501:600,1)); hold on
xlabel('Wrist Extension (°)')
ylabel('Active Wrist Extension Moment (Nm)')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')

%Passive extension moment
nexttile
plot(wrist_extension,br_moment_f(1:100,2)); hold on
plot(wrist_extension,br_moment_f(101:200,2)); hold on
plot(wrist_extension,br_moment_f(201:300,2)); hold on
plot(wrist_extension,br_moment_f(301:400,2)); hold on
plot(wrist_extension,br_moment_f(401:500,2)); hold on
plot(wrist_extension,br_moment_f(501:600,2)); hold on
xlabel('Wrist Extension (°)')
ylabel('Passive Wrist Extension Moment (Nm)')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')

%Active extension normalised
nexttile
plot(wrist_extension,br_moment_f(1:100,1)./br_moment_f(1:100,1)); hold on
plot(wrist_extension,br_moment_f(101:200,1)./br_moment_f(1:100,1)); hold on
plot(wrist_extension,br_moment_f(201:300,1)./br_moment_f(1:100,1)); hold on
plot(wrist_extension,br_moment_f(301:400,1)./br_moment_f(1:100,1)); hold on
plot(wrist_extension,br_moment_f(401:500,1)./br_moment_f(1:100,1)); hold on
plot(wrist_extension,br_moment_f(501:600,1)./br_moment_f(1:100,1)); hold on
xlabel('Wrist Extension (°)')
ylabel('Normalised Active Wrist Extension Moment')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')

%Passive extension normalised
nexttile
plot(wrist_extension,br_moment_f(1:100,2)./br_moment_f(1:100,2)); hold on
plot(wrist_extension,br_moment_f(101:200,2)./br_moment_f(1:100,2)); hold on
plot(wrist_extension,br_moment_f(201:300,2)./br_moment_f(1:100,2)); hold on
plot(wrist_extension,br_moment_f(301:400,2)./br_moment_f(1:100,2)); hold on
plot(wrist_extension,br_moment_f(401:500,2)./br_moment_f(1:100,2)); hold on
plot(wrist_extension,br_moment_f(501:600,2)./br_moment_f(1:100,2)); hold on
xlabel('Wrist Extension (°)')
ylabel('Normalised Passive Wrist Extension Moment')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')

%BR Wrist Deviation Moment at 30 degrees elbow flex
nexttile
plot(wrist_deviation,br_d(1:100,5)); hold on
plot(wrist_deviation,br_d(101:200,5)); hold on
plot(wrist_deviation,br_d(201:300,5)); hold on
plot(wrist_deviation,br_d(301:400,5)); hold on
plot(wrist_deviation,br_d(401:500,5)); hold on
plot(wrist_deviation,br_d(501:600,5)); hold on
plot(wrist_deviation,dev_m_pre);
xlabel('Radial Deviation (°)')
ylabel('Radial Deviation Moment (Nm)')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6','ECRB unparalysed')

%Normalised BR Wrist Deviation Moment at 30 degrees elbow flex
nexttile
plot(wrist_deviation,br_d(1:100,5)./dev_m_pre); hold on
plot(wrist_deviation,br_d(101:200,5)./dev_m_pre); hold on
plot(wrist_deviation,br_d(201:300,5)./dev_m_pre); hold on
plot(wrist_deviation,br_d(301:400,5)./dev_m_pre); hold on
plot(wrist_deviation,br_d(401:500,5)./dev_m_pre); hold on
plot(wrist_deviation,br_d(501:600,5)./dev_m_pre); hold on
plot(wrist_deviation,dev_m_pre./dev_m_pre)
xlabel('Radial Deviation (°)')
ylabel('Normalised Radial Deviation Moment')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6','ECRB unparalysed')


%---     ---    ---   Figure 3 (Active/ Passive force)  ---     ---     ---

figure('name','Active and Passive Force Components')
tiledlayout(2,3)

% Total Active/ Passive Split at 130 degrees elbow flex - SL1
nexttile
plot(wrist_extension,br_f_ap(1:100,10)); hold on
plot(wrist_extension,br_f_ap(1:100,11)); hold on
plot(wrist_extension,br_f_ap(1:100,12));
xlabel('Wrist Extension (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

% Total Active/ Passive Split at 130 degrees elbow flex - SL2
nexttile
plot(wrist_extension,br_f_ap(101:200,10)); hold on
plot(wrist_extension,br_f_ap(101:200,11)); hold on
plot(wrist_extension,br_f_ap(101:200,12));
xlabel('Wrist Extension (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

% Total Active/ Passive Split at 130 degrees elbow flex - SL3
nexttile
plot(wrist_extension,br_f_ap(201:300,10)); hold on
plot(wrist_extension,br_f_ap(201:300,11)); hold on
plot(wrist_extension,br_f_ap(201:300,12));
xlabel('Wrist Extension (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

% Total Active/ Passive Split at 130 degrees elbow flex - SL4
nexttile
plot(wrist_extension,br_f_ap(301:400,10)); hold on
plot(wrist_extension,br_f_ap(301:400,11)); hold on
plot(wrist_extension,br_f_ap(301:400,12));
xlabel('Wrist Extension (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

% Total Active/ Passive Split at 130 degrees elbow flex - SL5
nexttile
plot(wrist_extension,br_f_ap(401:500,10)); hold on
plot(wrist_extension,br_f_ap(401:500,11)); hold on
plot(wrist_extension,br_f_ap(401:500,12));
xlabel('Wrist Extension (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

% Total Active/ Passive Split at 130 degrees elbow flex - SL6
nexttile
plot(wrist_extension,br_f_ap(501:600,10)); hold on
plot(wrist_extension,br_f_ap(501:600,11)); hold on
plot(wrist_extension,br_f_ap(501:600,12));
xlabel('Wrist Extension (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

% Active Force at 130 degrees elbow flex 
nexttile
plot(wrist_extension,br_f_ap(1:100,10)); hold on
plot(wrist_extension,br_f_ap(101:200,10)); hold on
plot(wrist_extension,br_f_ap(201:300,10)); hold on
plot(wrist_extension,br_f_ap(301:400,10)); hold on
plot(wrist_extension,br_f_ap(401:500,10)); hold on
plot(wrist_extension,br_f_ap(501:600,10)); hold on
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')
xlabel('Wrist Extension (°)');
ylabel('Active Force (N)')

% Passive Force at 130 degrees elbow flex 
nexttile
plot(wrist_extension,br_f_ap(1:100,11)); hold on
plot(wrist_extension,br_f_ap(101:200,11)); hold on
plot(wrist_extension,br_f_ap(201:300,11)); hold on
plot(wrist_extension,br_f_ap(301:400,11)); hold on
plot(wrist_extension,br_f_ap(401:500,11)); hold on
plot(wrist_extension,br_f_ap(501:600,11)); hold on
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')
xlabel('Wrist Extension (°)');
ylabel('Passive Force (N)')

% Total Force at 130 degrees elbow flex 
nexttile
plot(wrist_extension,br_f_ap(1:100,12)); hold on
plot(wrist_extension,br_f_ap(101:200,12)); hold on
plot(wrist_extension,br_f_ap(201:300,12)); hold on
plot(wrist_extension,br_f_ap(301:400,12)); hold on
plot(wrist_extension,br_f_ap(401:500,12)); hold on
plot(wrist_extension,br_f_ap(501:600,12)); hold on
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')
xlabel('Wrist Extension (°)');
ylabel('Total Force (N)')

% --------------------------------------------------------------------------
% --------------------ELBOW FLEXION-----------------------------------------
% --------------------------------------------------------------------------

% -------------------Wrist Extension MOMENT-----------------------------------

figure('name','Effect of Elbow Position on Wrist Extension Moment')
tiledlayout(2,2)
nexttile()
plot(elbow_flexion,(br_el(1:100,2))); hold on
plot(elbow_flexion,(br_el(101:200,2))); hold on
plot(elbow_flexion,(br_el(201:300,2))); hold on
plot(elbow_flexion,(br_el(301:400,2))); hold on
plot(elbow_flexion,(br_el(401:500,2))); hold on
plot(elbow_flexion,(br_el(501:600,2))); hold on
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')
xlabel('Elbow Flexion (°)');
ylabel('Total Wrist Extension Moment (Nm)')

nexttile
plot(elbow_flexion,br_moment_el(1:100,2)); hold on
plot(elbow_flexion,br_moment_el(101:200,2)); hold on
plot(elbow_flexion,br_moment_el(201:300,2)); hold on
plot(elbow_flexion,br_moment_el(301:400,2)); hold on
plot(elbow_flexion,br_moment_el(401:500,2)); hold on 
plot(elbow_flexion,br_moment_el(501:600,2)); hold on 
xlabel('Elbow Flexion (°)')
ylabel('Passive Wrist Extension Moment (Nm)')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')

nexttile
plot(elbow_flexion,br_moment_el(1:100,1)); hold on
plot(elbow_flexion,br_moment_el(101:200,1)); hold on
plot(elbow_flexion,br_moment_el(201:300,1)); hold on
plot(elbow_flexion,br_moment_el(301:400,1)); hold on
plot(elbow_flexion,br_moment_el(401:500,1)); hold on
plot(elbow_flexion,br_moment_el(501:600,1)); hold on
xlabel('Elbow Flexion (°)')
ylabel('Active Wrist Extension Moment (Nm)')
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6')

nexttile
plot(elbow_flexion,(br_el(1:100,2))); 
title('Slack Length 1');
xlabel('Elbow Flexion (°)');
ylabel('Wrist Extension Moment (Nm)') 

nexttile
plot(elbow_flexion,(br_el(101:200,2))); 
xlabel('Elbow Flexion (°)');
ylabel('Wrist Extension Moment (Nm)');
title('Slack Length 2');

nexttile
plot(elbow_flexion,(br_el(201:300,2))); 
xlabel('Elbow Flexion (°)');
ylabel('Wrist Extension Moment (Nm)')
title('Slack Length 3');

nexttile
plot(elbow_flexion,(br_el(301:400,2))); 
xlabel('Elbow Flexion (°)');
ylabel('Wrist Extension Moment (Nm)');
title('Slack Length 4');

nexttile
plot(elbow_flexion,(br_el(401:500,2))); 
xlabel('Elbow Flexion (°)');
ylabel('Wrist Extension Moment (Nm)');
title('Slack Length 5');

%--------------------ELBOW EXTENSION MOMENT--------------------------------

% Elbow Extension Moment

figure('name','Effect of Elbow Position on Elbow Flexion Moment')
tiledlayout(2,2)

nexttile()
plot(elbow_flexion,(br_el(1:100,3))); hold on
plot(elbow_flexion,(br_el(101:200,3))); hold on
plot(elbow_flexion,(br_el(201:300,3))); hold on
plot(elbow_flexion,(br_el(301:400,3))); hold on
plot(elbow_flexion,(br_el(401:500,3))); hold on
plot(elbow_flexion,(br_el(501:600,3))); hold on
plot(elbow_flexion,br_m_pre)
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6','Pre-transfer')
xlabel('Elbow Flexion (°)');
ylabel('Br Elbow Flexion Moment (Nm)')

nexttile()
plot(elbow_flexion,el_flex(1:100)); hold on
plot(elbow_flexion,el_flex(101:200)); hold on
plot(elbow_flexion,el_flex(201:300)); hold on
plot(elbow_flexion,el_flex(301:400)); hold on
plot(elbow_flexion,el_flex(401:500)); hold on
plot(elbow_flexion,el_flex(501:600)); hold on
plot(elbow_flexion,el_m_pre)
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6','Pre-transfer')
xlabel('Elbow Flexion (°)');
ylabel('Total Elbow Flexion Moment (Nm)')

nexttile()
plot(elbow_flexion,el_flex(1:100)./el_m_pre); hold on
plot(elbow_flexion,el_flex(101:200)./el_m_pre); hold on
plot(elbow_flexion,el_flex(201:300)./el_m_pre); hold on
plot(elbow_flexion,el_flex(301:400)./el_m_pre); hold on
plot(elbow_flexion,el_flex(401:500)./el_m_pre); hold on
plot(elbow_flexion,el_flex(401:500)./el_m_pre); hold on
plot(elbow_flexion,el_m_pre./el_m_pre)
legend('Length 1','Length 2','Length 3','Length 4','Length 5','Length 6','Pre-transfer')
xlabel('Elbow Flexion (°)');
ylabel('Normalised Total Elbow Flexion Moment (Nm)')

nexttile
plot(elbow_flexion,(br_el(1:100,3)));
title('Slack length 1')
xlabel('Elbow Flexion (°)');
ylabel('Elbow Flexion Moment (Nm)')

nexttile
plot(elbow_flexion,(br_el(101:200,3))); 
title('Slack length 2')
xlabel('Elbow Flexion (°)');
ylabel('Elbow Flexion Moment (Nm)')

nexttile
plot(elbow_flexion,(br_el(201:300,3)));
title('Slack length 3')
xlabel('Elbow Flexion (°)');
ylabel('Elbow Flexion Moment (Nm)')

nexttile
plot(elbow_flexion,(br_el(301:400,3))); 
title('Slack length 4')
xlabel('Elbow Flexion (°)');
ylabel('Elbow Flexion Moment (Nm)')

nexttile
plot(elbow_flexion,(br_el(401:500,3)));
title('Slack length 5')
xlabel('Elbow Flexion (°)');
ylabel('Elbow Flexion Moment (Nm)')


%-----------------ACTIVE/ PASSIVE FORCE SPLIT------------------------------

% Active/ Passive Force
figure('name','Effect of Elbow Position on Active and Passive Force')
nexttile
plot(elbow_flexion,br_el_ap(1:100,1)); hold on
plot(elbow_flexion,br_el_ap(1:100,2)); hold on
plot(elbow_flexion,br_el_ap(1:100,3));
title('Brachioradialis Force - Slack Length 1')
xlabel('Elbow Flexion (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

nexttile
plot(elbow_flexion,br_el_ap(101:200,1)); hold on
plot(elbow_flexion,br_el_ap(101:200,2)); hold on
plot(elbow_flexion,br_el_ap(101:200,3));
title('Brachioradialis Force - Slack Length 2')
xlabel('Elbow Flexion (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

nexttile
plot(elbow_flexion,br_el_ap(201:300,1)); hold on
plot(elbow_flexion,br_el_ap(201:300,2)); hold on
plot(elbow_flexion,br_el_ap(201:300,3));
title('Brachioradialis Force - Slack Length 3')
xlabel('Elbow Flexion (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

nexttile
plot(elbow_flexion,br_el_ap(301:400,1)); hold on
plot(elbow_flexion,br_el_ap(301:400,2)); hold on
plot(elbow_flexion,br_el_ap(301:400,3));
title('Brachioradialis Force - Slack Length 4')
xlabel('Elbow Flexion (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

nexttile
plot(elbow_flexion,br_el_ap(401:500,1)); hold on
plot(elbow_flexion,br_el_ap(401:500,2)); hold on
plot(elbow_flexion,br_el_ap(401:500,3));
title('Brachioradialis Force - Slack Length 5')
xlabel('Elbow Flexion (°)')
ylabel('Force (N)')
legend('Active','Passive','Total')

nexttile
plot(elbow_flexion,br_el_ap(1:100,1)); hold on
plot(elbow_flexion,br_el_ap(101:200,1)); hold on
plot(elbow_flexion,br_el_ap(201:300,1)); hold on
plot(elbow_flexion,br_el_ap(301:400,1)); hold on
plot(elbow_flexion,br_el_ap(401:500,1)); hold on
title('Active Brachioradialis Force')
legend('Length 5','Length 6','Length 7','Length 8','Length 9')
xlabel('Elbow Flexion (°)');
ylabel('Active Force (N)')

nexttile
plot(elbow_flexion,br_el_ap(1:100,2)); hold on
plot(elbow_flexion,br_el_ap(101:200,2)); hold on
plot(elbow_flexion,br_el_ap(201:300,2)); hold on
plot(elbow_flexion,br_el_ap(301:400,2)); hold on
plot(elbow_flexion,br_el_ap(401:500,2)); hold on
title('Passive Brachioradialis Force')
legend('Length 5','Length 6','Length 7','Length 8','Length 9')
xlabel('Elbow Flexion (°)');
ylabel('Passive Force (N)')

nexttile
plot(elbow_flexion,br_el_ap(1:100,3)); hold on
plot(elbow_flexion,br_el_ap(101:200,3)); hold on
plot(elbow_flexion,br_el_ap(201:300,3)); hold on
plot(elbow_flexion,br_el_ap(301:400,3)); hold on
plot(elbow_flexion,br_el_ap(401:500,3)); hold on
title('Total Brachioradialis Force')
legend('Length 5','Length 6','Length 7','Length 8','Length 9')
xlabel('Elbow Flexion (°)');
ylabel('Total Force (N)')


% Total Moment
figure('name','total moment')
tiledlayout('flow')
plot(wrist_extension,total_moment(1:100)); hold on
plot(wrist_extension,total_moment(101:200)); hold on
plot(wrist_extension,total_moment(201:300)); hold on
plot(wrist_extension,total_moment(301:400)); hold on
plot(wrist_extension,total_moment(401:500)); hold on
plot(wrist_extension,total_m_pre)
legend('1','2','3','4','5','6')

nexttile
plot(wrist_extension,br_f(1:100,5)); hold on
plot(wrist_extension,br_f(101:200,5)); hold on
plot(wrist_extension,br_f(201:300,5)); hold on
plot(wrist_extension,br_f(301:400,5)); hold on
plot(wrist_extension,br_f(401:500,5)); hold on
plot(wrist_extension,ecrb_m_pre);
title('Total Wrist Extension Moment','30° Elbow Flexion')
xlabel('Wrist Extension (°)')
ylabel('Wrist Extension Moment (Nm)')
legend('Length 5','Length 6','Length 7','Length 8','Length 9','pre')



