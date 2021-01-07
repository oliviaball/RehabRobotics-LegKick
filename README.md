# RehabRobotics-LegKick

The goal of the assignment was to develop a simulation of a leg exoskeleton for a hip, knee and ankle in the sagittal plane. The follow tasks were used to take the given information and turn it into joint torques that were required to produce that motion.

1) Create a model of each limb and device to determine the location of the combined center of mass and moment of inertia. The exoskeleton had forces of 50 N, 25 N and 25 N located at the proximal ends of the hip, knee and ankle respectively. 

	•	Open the file desiredtorques.m.
	•	Ensure the file leg_kick_rad.xlsx is also downloaded and in the same folder. This was a given document that contains the desired angle, velocity, and acceleration for each of the three segments.
	•	Take own measurements for length and radius of each segment as described and fill them in. This ends on line 22 and no other modifications need to be made.
	•	The command window will tell you the location of the combined center of mass and the moment of inertia for each segment. 

2) Using given file leg_kick.xls, compute the desired torque needed at each joint at each time step to create the desired torque at each joint to create the desired movement. The file will be saved in a text file where the first column corresponds to joint 1 (upper leg), joint 2 (lower leg), and joint 3 (ankle). 

	•	The text file will be saved as lastname_torquevals.txt - these values will be used to control the motors in Working Model.

3) Create a text file which includes the desired joint angles.  

	•	Open create_traj.m - the file will create a trajoutput.txt file which is the desired trajectory based on the length of your leg. DO NOT CLEAR THE WORKBOOK - it will automatically include the values for your leg length.
	•	Next run the file desiredjointangles.m which will create a final text file - lastname_jointsngles.txt

4) The text files were then used to create the video that is included. The model works as expected. 

5) The numbers within desiredtorques.m can be modified with your leg measurements to ensure the model works.
