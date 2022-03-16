# human-walking-biomechanics

## Summary ##

Force plate data and motion capture data were collected for 33 different walking trials with 
various step lengths, step widths, step frequencies, and speeds. Ground reaction force data were 
sampled at 1200 Hz and motion capture data were sampled at 120 Hz. The data files contain variables that were
created and exported from Visual 3D. For more information about all the variables see Variable 
Descriptions below.

### Data Files ###

There are two versions of the data:

- The all strides data files is all the data exported from V3D (the entire trial)

- The 5 strides data files contains data for 5 selected strides

There are 10 raw data files, one for each of 10 subjects. For example, for subject 1, the raw data filename is 
**p1_AllStridesData.mat**. This file contains a struct named data containing 9 fields. The variables contained within the
fields are described below in Variable Descriptions. 

There are 9 5-step files since there is no 5-step file for subject 10 (see Missing/Excluded Data). For subject 2 the 5 step 
data file name is **p2_5StridesData.mat**. These files contain a struct called data with the same 9 fields as the raw data. 
All of the variables are the same as the raw data except instead of containing data for the entire trial it only has 
data for 5 strides. The interval for the 5 strides was selected by looking at ground reaction force plots and joint power
plots. Intervals without any missteps on the force platform were chosen.

### Trials ###

There are a total of 33 trials that can be divided into the following 5 experimental conditions:

1. Constant step length (increasing speed),

2. Constant step frequency (increasing speed)

3. Constant speed (increasing step length)

4. Preferred walking (increasing speed)

5. Constant speed (increasing step width)

See Trial Descriptions below for more detail.


## Variable Descriptions ##

Here are some abbreviations used in the variable names:

|  Variable   | Full Names/Description |
| ----------- |		-----------		       |
|    l/L      |    Left 	             |
|    r/R      |    Right    	         |
| ----------- | ---------------------- |
|     ft      |    Foot	          
|     sk      |    Shank   	       
|     th      |    Thigh	             |
|     pv      |    Pelvis    	         |
| ----------- | ---------------------- |
|    ank      |    ankle 	             |
|    kne      |    knee	      		     |
|    hip      |    hip 	               |
| ----------- | ---------------------- |
|    loc      | local coordinates      |
|    lab      |  lab coordinates       |
|   dist      |    distal    	         |
|    prox     |    proximal    	       |
| ----------- | ---------------------- |
|    wrt	    | with respect to	       |
|    cop	    | centre of pressure     |
| ----------- | ---------------------- |
|   original  |    Original Data  	   |
|processed/proc|    Processed Data     |


### **Force Platform Variables** ###

*Location*: These variables are stored in field called ‘Platform’

*Variables*:

•	**ForcePlatformOrigin**- float variable containing force platform origin in the force platform coordinate system, origin is the centre of the plate

•	**ForcePlatformCorners**- Lab coordinates defining force platform corners, the order is (+x, +y), (-x, +y), (-x, -y), (+x, -y)

•	**ForcePlatformCalibration**- Used to convert analog signals from Volts into Newtons 

•	**ForcePlatformZero**- frames that background noise is calculated, this is then subtracted from each analog channel

•	**ForcePlatformZeros**- Force platform zeroing range

•	**MatrixStore**- ForcePlatformCalibration array storage (variable= BYROW)

Example: **Data.ForcePlatform.ForcePlatformZero**

<https://www.c-motion.com/v3dwiki/index.php?title=C3D_Format#FORCE_PLATFORM_Parameters>


### **Analog** ###

*Location*: These variables are stored in field called ‘Analog’

*Variables*: These variables start with an **f** or **m** followed by a **1** or **2**, and finally are suffixed with **original** or **processed**. 
These variables contains raw data from one of the 2 force plates, in all 3 directions (columns 1, 2, and 3 are the x, y, and z directions, respectively). 

*Description*:

•	**f**- ground reaction force

•	**m**- moment

•	**1**-	Left force plate 

•	**2-**	Right force plate

•	**original**- not processed ex. f1zoriginal

•	**processed**- processed ex. m2yprocessed

*Coordinate system*: lab

*Units*: Volts

*Processing*: The processed data is processed using a butterworth filter and a frequency cut-off of 25.

Example: **Data.Analog.f1original**


### **Target Data** ###

*Location*: These variables are stored in field called ‘TargetData’

*Variables*: Variables are prefix with **l** for left and **r** for right, followed by the anatomical landmark, and suffixed with **_pos** 
for the raw data and **_pos_proc** for the processed data. 

*Anatomical Landmark*: 

•	**5TH**- marker on 5th toe  

•	**AC**- acromion process

•	**ASI**- anterior superior iliac crest

•	**CAL**- calcaneus

•	**EP**- elbow epicondyle 

•	**GTR**- greater trochanter

•	**LEP**- lateral knee epicondyle

•	**LML**- lateral malleolus (outside of ankle)

•	**MEP**- Medial knee epicondyle

•	**MML**- medial malleolus (inside of ankle)

•	**SH1**- lateral shank marker 1

•	**SH2**- lateral skank marker 2

•	**SH3**- lateral shank marker 3

•	**TH1**- lateral thigh marker 1

•	**TH2**- lateral thigh marker 2

•	**TH3**- lateral thigh marker 3

•	**WR**- wrist

•	**SACR**-sacrum (not prefixed with L or R as this marker is centred)

*Orientation*: Lab coordinates

*Units*: metres

*Processing*: Processed with a butterworth filter, Num_Reflected value of 120 and a frequency cut-off of 6. Sampled at 120 Hz.

*Description*: The variables have five columns: The first three are the x, y, z, positions in space. The fourth is the component 
with info from the motion capture hardware regarding reliability. The fifth column is the camera contribution, it is telling us how 
many cameras were used for the reconstruction. The processed variables have zeros for the fourth and fifth column so look to the 
original variables for this information.

Example: **Data.TargetData.L5TH_pos** contains raw data and L5TH_pos_proc contains processed data for the marker on the 5th toe.


### **Landmarks** ###

*Location*: These variables are stored in field called ‘Landmark’

*Description*: 	**HHL** and **HHR** are the left and right hip joints centres. The location of the hip joint centre is calculated using the pelvis markers and 
the Helen Hayes pelvis segment model.

*Units*: metres (m)

<https://www.c-motion.com/v3dwiki/index.php/Helen_Hayes_(Davis)_Pelvis>


### **Time** ###

*Location*: These variables are stored in field called ‘Time’

*Variables*:

•	**AnalogTime**- Time recorded at 1200 Hz (force platform data)

•	**FRAMES**- frame number

•	**TIME**- Time recorded at 120 Hz (motion capture data)

Example: **Data.Time.FRAMES**


### **Ground Reaction Force Data** ###

*Location*: These variables are stored in field called ‘Force’

*Coordinate System*: lab

*Variables*: **force1**, **force2**, **cop1**, **cop2**, **freemoment1**, and **freemoment2**

**1**-	Left force plate      **2**-	Right force plate

•	**force**- force data for the x, y, and z directions is stored in the 1st, 2nd, and 3rd columns

*Units*: Newtons (N)

*Calculated*:  by subtracting the background noise (stored in ForcePlatformZero variable) from the analog signal and then multiplying the analog force variables
by the calibration matrix (ForcePlatformCalibration). This converts the signal in volts into a force with units Newton’s and removes background noise.

•	**cop**- centre of pressure data for the x, y, and z directions is stored in the 1st, 2nd, and 3rd column

*units*: metres (m)

*Calculated*: using the ForcePlatformOrigin variable and the forec1 and force2 variables explained above. It also uses the analog moment variables (M) 
multiplied by the calibration matrix (ForcePlatformCalibration). The COP calculations are: 

COP[X]= (ForcePlatformOrigin[Z] * Force[X] – My) / Force[Z]

COP[Y]= (ForcePlatformOrigin[Z] * Force[Y] – Mx) / Force[Z]

COP[Z]= ForcePlatformOrigin[Z]

•	**freemoment**- moment caused about the lab z-axis due to the ground reaction force

*Calculated*: using the cop variables above as well as the analog force and moment varibales multiplied by the calibration matix (ForcePlatformCalibration). 
The freemoment calculations are:

freemomnt[X]= 0

freemoment[Y]= 0

freemoment[Z]= Mz – (COP[X] * Force[Y] – COP[Y] * Force[X])

Example: **Data.Force.cop2** contains the centre of pressure for the right side.


### **Kinetic_Kinematic Data** ###

*Location*: These variables are stored in field called ‘Kinetic_Kinematic’

*Notes*: These 12 variables can have the prefix **lFt, lSh, lTh, rFt, rSh, rTH, or rPV**.  Ex. **lFtAngAcc**: left foot angular acceleration

*Variable Types*: 

•	**AngAcc**- Angular acceleration of the segment (rad/s2)

•	**AngVel**- Angular Velocity of the segment (rad/s)

•	**CGAcc**- Segment Center of Mass Acceleration (m/s2)

•	**CGPos** - Segment Center of Mass Position (m)

•	**CGVel** - Segment Center of Mass Velocity (m/s)

•	**DistEndPos**- Position of the Distal End of the segment (m)

•	**DistEndVel** – Velocity of the Distal End of the segment (m/s)

•	**ProxEndForce** - Proximal End Joint Reaction Force (N).  This is the force acting on the proximal end of the joint, for example, lthProxEndForce is the force acting from the hip on the proximal end of the thigh and thus the z-component of the joint reaction force will be primarily negative values

•	**ProxEndPos** - Proximal End Position (m)

•	**ProxEndTorque** - Proximal End Joint Moment (Nm)

•	**ProxEndVel** - Proximal End Velocity (m/s)

•	**SegResidual** - the result of the least squares fit of the tracking markers at each frame of data

All variables have 3 columns of data except for the SegResidual variables (they have 1 column of data)
 
*Coordinate System*: Lab

Example: **Data.Kinetic_Kinematic.rShDistEndVel**

<https://www.c-motion.com/v3dwiki/index.php?title=KINETIC_KINEMATIC>


### **Link_Model_Based Variables** ###

*Location*: These variables are stored in field called ‘Link_Model_Based’

*Notes*: 

1.	These variables all have 3 columns of data acting along the x, y, and z, directions or for data rotating about these 3 axes. 

2.	The joint force, joint moment, and joint power variables, they refer to the proximal segment ex. right ankle force is resolved in the right shank coordinate system

3. Shortened joint names: ank, kne, and hip.

•	**Joint Angle**

*Description*: 	Joint Angle variables follow the pattern of **l** or **r** for left or right, the shortened joint name, and suffix with **angle**.  
Variable defines the difference in rotation of the proximal and distal segment in the proximal segments coordinate system. These variables are normalized 
so that 0° occurs when the subject is in the reference position. 

*Example*: **l_ank_angle** for left ankle angle.

For the Ankle: dorsiflexion, ffadduction, and eversion are positive- adduction/abduction is also referred to as internal/external rotation

For the Knee: extension, adduction and internal rotation are positive

For the Hip: flexion, adduction, and internal rotation are positive

*Units*: degrees

•	**Joint Moment**

*Description*: 	Joint Moment variables follow the pattern of **l** or **r** for left or right, the shortened joint name, and suffix with **moment**. 

*Example*: **r_kne_moment** for right knee moment

For the Ankle: moments acting in the direction of dorsiflexion, ffadduction, and eversion are positive- adduction/abduction is also referred to as 
internal/external rotation

For the Knee: moments causing extension, adduction and internal rotation are positive

For the Hip: moments in the direction of flexion, adduction, and internal rotation are positive

*Notes*: These variables are normalized using subjects’ body mass.

*Units*: Nm/kg

•	**Joint Power**

*Description*: 	Joint Power variables follow the pattern of **l** or **r** for left or right, the shortened joint name, and suffix with **power**. 

*Example*: **r_hip_power** for right hip power

This variable has 3 column of scaler data since it is the dot product of the joint moment and the angular velocity resolved in the proximal segments 
coordinate system

*Notes*: These variables are normalized using the subjects’ body mass

*Units*: W/kg

•	**Joint Force**

*Description*: 	Joint Force variables follow the pattern of **l** or **r** for left or right, the shortened joint name, and suffix with **force_loc**. 
These variables are the joint reaction force at the proximal end of a segment acting from the segment to the joint

*Example*: **r_hip_force_loc** for right hip force. 

*Notes*: Normalized using subject body mass

*Units*: N/kg

•	**Joint Velocity**

*Description*: Joint Velocity variables follow the pattern of **l** or **r** for left or right, the shortened joint name, and suffix with **vel**. 
This variable describes the angular velocity of the distal segment relative to the proximal segment. 

*Example*: **r_ank_vel** for right ankle velocity.

*Units*: degrees/s 

•	**Segment RotEnergy**

*Description*: Segment RotEnergy variables follow the pattern of **l** or **r**, followed by the shorten segment name, and finally, suffix with **_rotenergy**. 

*Example*: **rft_rotenergy** is the rotational energy of the right foot.

*Calculation*: Rotational Energy= 0.5*IXX*AngVel.X*AngVel.X + 0.5*IYY*AngVel.Y*AngVel.Y + 0.5*IZZ*AngVel.Z*AngVel.Z

*Units*: Joules (J)

•	**Segment Position wrt Distal Joint**

*Description*: 	Segment Position variables follow the pattern **l** or **r**, followed by shorten segment name and **wrt**, then **l** or **r**, followed by 
the reference segment’s shorten name. 

*Example*: **lsk_wrt_lank** is the position of the left shank with reference to the left ankle.

*Coordinate System*: 	The proximal segment centre of gravity is in the distal segment’s coordinate system

*Units*: metres


### **Subject Variables** ###

*Location*: These variables are stored in the field called ‘Subject’

*Variables*: 

•	**Height**- Subject’s height in metres

•	**Mass**- Subject’s mass in kilograms

•	**centerofmass**- Subject’s centre of mass in metres

Example: **Data.Subject.Mass** contains the mass of the subject.


## Trial Descriptions ##

Trials correspond to the following conditions:

1.	0.7 m/s constant step length
2.	0.7 m/s preferred
3.	0.7 m/s constant step frequency
4.	0.9 m/s constant step length
5.	0.9 m/s preferred
6.	0.9 m/s constant step frequency
7.	1.1 m/s constant step length
8.	1.1 m/s preferred
9.	1.1 m/s constant step frequency
10.	1.6 m/s constant step frequency
11.	1.6 m/s preferred
12.	1.6 m/s constant step length
13.	1.8 m/s constant step frequency
14.	1.8 m/s preferred
15.	1.8 m/s constant step length
16.	2.0 m/s preferred
17.	1.25 m/s lowest step frequency
18.	1.25 m/s lower step frequency
19.	1.25 m/s low step frequency
20.	1.25 m/s preferred
21.	1.25 m/s preferred
22.	1.25 m/s preferred
23.	1.25 m/s high step frequency
24.	1.25 m/s higher step frequency
25.	1.25 m/s highest step frequency
26.	1.25 m/s zero step width
27.	1.25 m/s 10 cm step width
28.	1.25 m/s 20 cm step width
29.	1.25 m/s 30 cm step width
30.	1.25 m/s 40 cm step width
31.	1.40 m/s constant step length
32.	1.40 m/s preferred
33.	1.40 m/s constant step frequency

In the above:

•	preferred means using preferred step length & step frequency for the speed in question

•	constant step length means using the step length that was preferred when walking at 1.25 m/s

•	constant step frequency means using the step frequency that was preferred when walking at 1.25 m/s


## Missing/Excluded Data ##

### **All Strides Data** ###
1.	Subject 6 trial 31, the raw data will have variables full of zeros for the motion capture variables, and empty variables for all variables relying on motion capture data as there was no motion capture data collected for this trial.
2.	Subject 7 trial 24 doesn’t exist
3.	Subject 9 does not have upper body motion capture data for any of the trials. The following markers were not used for Subject 9: LAC, RAC, LEP, REP, LWR, RWR. Any variables related to these variables are empty in the raw data.

### **5 Step Data** ###
1.	For all subjects, trials 26-30 are not included in the 5-step data. These trials are the step width trials and since the heel strikes are difficult to identify, defining steps and picking a 5-step range was not complete.
2.	Subject 3 trials 4 has a sync error. It is included in the raw data but not in the 5-step data.
3.	Subject 4 trial 1 there was incorrect stepping on the force platform so 5 step range was not selected.
4.	Subject 6 is missing 5 step data for trial 21 and 31 as no 5-step range could be selected. For trial 21 there is only data for 4 steps and for trial 31 there is no motion capture data. For trial 31, the raw data will have variables full of zeros for the motion capture variables, and empty variables for all variables relying on motion capture data.
5.	Subject 7 trials 24 doesn’t exist.
6.	Subject 9 does not have upper body motion capture data for any of the trials. The following markers were not used for Subject 9: LAC, RAC, LEP, REP, LWR, RWR. Any variables related to these variables do not exist in the 5-step data.
7.	Subject 9 trials 14 has a sync error. It is included in the raw data but not in the 5-step data.
8.	No 5-step data for Subject 10. Subject 10 has small pelvis measurements, this is possibly a result of incorrect marker placement. Due to this error, 5-step intervals could not be accurately selected for this subject.




