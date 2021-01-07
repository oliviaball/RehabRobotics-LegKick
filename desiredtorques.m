% Pull in angle, velocioty and acceleration for desired movement 
legkick= xlsread('leg_kick_rad.xlsx'); %CHANGE IF NAME IS DIFFERENT
t1=legkick(:,2); % angle for the upper leg (rad)
td1=legkick(:,3); % velocity for the upper leg (rad/sec)
tdd1=legkick(:,4); % acceleration for the upper leg (rad/sec^2)
t2=legkick(:,5); % angle for the lower leg (rad)
td2=legkick(:,6); % velocity for the lower leg (rad/sec)
tdd2=legkick(:,7); % acceleration for the lower leg (rad/sec^2)
t3=legkick(:,8); % angle for the foot (rad)
td3=legkick(:,9); % velocity for the upper leg (rad/sec)
tdd3=legkick(:,10); % acceleration for the foot(rad/sec^2)
t4=zeros(1,length(t1))'; %angle doesnt change between foot and ankle so set = 0 
torquevals=[zeros(1,101)' zeros(1,101)' zeros(1,101)']; %for indexing all the torques at the end 

%*******************FILL OUT THESE FOR OWN LEG *******************% 
%generalized enough that you can adjust these and still run the model for people of different heights  
lu=0.42; %m - (lu = length of upper leg), (16.5 inches)
ll=0.38; %m - (ll = length of lower leg), (15 inches)
lf=0.254; %m - (lf = length foot), (10 inches)
ru=0.089; %m - (ru = radius upper), D=7 inches, R=3.5 inches
rm=0.0635; %m - (rm = radius middle or radius of the knee), D=5 inches, R=2.5 inches
rl=0.04445; %m - (rl = radius lower leg or ankle area), D=3.5 inches, R=1.75 inches
rf=0.035; %m - (rf = radius foot), D=2.75 inches, R=1.375 inches

% Pregiven assumtions - NO ACTION NEEDED
lt_x=rf/2; %m - (lt_x = half length of foot from heel to toes)
lt_y=lf/2; %m - (lt_y = half of length from top to bottom of foot)
d=1100; %kg/m^3, assumed value for density
nmu=5; %kg - point force 50 N 
nml=2.5; %kg - point force 25 N
nmf=2.5; %kg - point for 25 N

%************** Mass and Center of Mass *************
% Mass and Location New COM - UPPER LEG
[mu, cmu, mcmu] = cylcoms (ru,rm,lu,d,nmu);
% Mass and Location New COM - LOWER LEG
[ml,cml, mcml] = cylcoms (rm,rl,ll,d,nml);
% Mass and Location New COM - FOOT
[mf,cmf, mcmf] = cylcoms (rm,rl,ll,d,nml);
%******ONLY FOR THE FOOT SINCE HAS Y-DIRECTION TO ACCOUNT FOR ********
cmf_x=((mf*lt_x+nmf*0)/(mf+nmf)); %distance of combined COM in m
cmf_y=((mf*(lt_y-rl)+nmf*0)/(mf+nmf)); %distance in the y direction for COM, use (lt_y-rl) to get distance from center of the foot to the point mass that is in the center of the foot 

% Print out the values 
fprintf('The location of the combined center of mass from the proximal end of the upper leg is:%d.\n',cmu)
fprintf('The location of the combined center of mass from the proximal end of the lower leg is:%d.\n',cml)
fprintf('The location of the combined center of mass from the proximal end of the foot is:%d.\n',cmf)

%************** MOMENT OF INERTIA *************
% Moment of Inertia - UPPER LEG
[ixxu, iyyu, izzu] = momentinter(ru,rm,cmu,mcmu);
% Moment of Inertia - LOWER LEG
[ixxl, iyyl, izzl] = momentinter(rm,rl,cml,mcml);
% Moment of Inertia - FOOT
[ixxf, iyyf, izzf] = momentinter(rl,rf,cmf,mcmf);

% Print out the values 
fprintf('The moment of inertia in the z plane for the upper leg is:%d.\n',izzu)
fprintf('The moment of inertia in the z plane for the lower leg is:%d.\n',izzl)
fprintf('The moment of inertia in the z plane for the foot is:%d.\n',izzf)

%% Finding the desired torque at each input 
for i=1:101
%Collect the inputs, the for loop will loop through each of the timesteps
ta=[t1 t2 t3]; % Collect all of the angles (radians)
thetadota=[td1 td2 td3]; % Collect all of the velocities (radians/second^2)
thetadda= [tdd1 tdd2 tdd3]; % Collect all of the accelerations rad^2/sec)

% Pull theta 1, theta 2 and theta 3 at each timestep 
t=ta(i,:);
thetadot=thetadota(i,:); 
thetadd=thetadda(i,:);

% Outward recursion important conditions:
p=[0 0 0;lu 0 0;ll 0 0]'; %position from transform matrix 
s=[cmu cml cmf_x; 0 0 cmf_y; 0 0 0 ]; % Scaling conditions
m=[mcmu mcml mcmf]; % Mass conditions in kilograms
ic=[0 0 0; 0 0 0; izzu izzl izzf]; % Inertia tensor stress

%Iinward recursion important conditions: 
p2=[lu 0 0;ll 0 0;lt_x lt_y 0]'; %iniital vector to include the gravity
s2=[cmu cml cmf_x; 0 0 cmf_y; 0 0 0]; %initial s conditions
x2a=[t2 t3 t4]; % Theta2 Theta 3 Thetainitial (Theta Initial assumed 0)
x2=x2a(i,:);

[N,Farray]=outwardthree(t,thetadot,thetadd,p,s,m,ic);
torque=inward(N,Farray,p2,s2,x2);

torquevals(i,:)=torque; % Gather all the desired data 
end 

% *************** Create a text file & plot ************************ %
time=0:0.01:1; %1 second, create for plotting
timeandtorque = [time' torquevals];

txtfile = fopen('lastname_torquevals.txt','w');
fprintf(txtfile,'%5f %5f %5f %5f\n',[timeandtorque].');
fclose(txtfile);


% torque1=torquevals(:,1);
% torque2=torquevals(:,2);
% torque3=torquevals(:,3);
% figure(1)
% plot(time,torque1)
% hold on
% plot(time,torque2)
% hold on
% plot(time,torque3)
% grid on
% legend('Joint 1', 'Joint 2', 'Joint 3')
% title('Plot of the Desired Torques for Each Joint')
% xlabel('Time (Seconds)')
% ylabel('Magnitude (N*m)')
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,jointforce]=outwardthree(t,thetadot,thetadd,p,s,m,ic)
%************** MATLAB "outwardthree function (Olivia Ball) *************
% The function completes the outward recursion for computing static
% endpoints forces and moments (Newton-Euler Inverse Dynamics Algorithm)

% SYNTAX        [N,Farray]=outwardthree(t,thetadot,thetadd,p,s,m,ic)
% INPUTS:       t          Angle desired  
%               thetadot   Velocity Desired 
%               thetadd	   Acceleration Desired
%               p          Position from Transform Matrix
%               s          Scaling Conditions 
%               m          Mass 
%               ic         Intitial Tensor Stress
%             
% OUTPUT:      N                  
%              jointforce         Joint Force  
%              
% CALLS:        
% CALLED BY:	
%~~~~~~~~~~~~~~~~~~~~~~ Begin Program: ~~~~~~~~~~~~~~~~~~~~~~~~~~
for j=2:4
    
    %Set the initial conditions
    g= -9.81; %m/sec^2m - gravity pulling down
    w(:,1)=[0;0;0]'; % Set inital condition for no rotation
    wd(:,1)=[0;0;0]'; % Set initial conditions for no accelerastion
    vd(:,1)=[g;0;0]'; %set init656ial condition for velocity 
    
    %Euler-Outward Recursion
    r=[cos(t(j-1)) sin(t(j-1)) 0;-sin(t(j-1)) cos(t(j-1)) 0; 0 0 1]; % Matrix for Ri,i-1
    wnew=r*w(:,j-1)+[0 0 thetadot(j-1)]';
    w(:,j)=[wnew]; 
    wdnew=r*wd(:,j-1)+cross(r*w(:,j-1),[0 0 thetadot(:,j-1)]')+[0 0 thetadd(:,j-1)]';
    wd(:,j)=[wdnew];
    cwp=cross((w(:,j-1)),(p(:,j-1)));
    multip=cross(wd(:,j-1),p(:,j-1))+cross(w(:,j-1),cwp);
    vdnew=r*((vd(:,j-1))+multip);
    vd(:,j)=[vdnew]; 
    cws=cross((w(:,j)),(s(:,j-1))); 
    vcdot= (vd(:,j))+(cross(wd(:,j),s(:,j-1)))+(cross(w(:,j),cws));
    vc(:,j)=[vcdot];
    F=m(:,j-1)*vc(:,j);
    jointforce(:,j)=[F];
    icw=cross(ic(:,j-1),w(:,j));
    n=ic(:,j-1).*wd(:,j)+cross(w(:,j),icw);
    N(:,j)=[n];
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function torque=inward(N,jointforce,pbase,sbase,xbase)
%************** MATLAB "mementinter function (Olivia Ball) *************
% Use of inward recursion to find the desired torque

% SYNTAX        torque=inward(N,Farray,p2,s2,x2)
% INPUTS:       N           Radius of one of the end point 
%               jointforce  Joint Force 
%               pbase       Position from the base
%               sbase       Scaling conditions 
%               xbase       Desired location
%             
% OUTPUT:      torque      Desired torque
%              
% CALLS:       
% CALLED BY:	
%~~~~~~~~~~~~~~~~~~~~~~ Begin Program: ~~~~~~~~~~~~~~~~~~~~~~~~~~

 for i=4:-1:2 
    
    % Set the initial conditions
    fin(:,4)=[0 0 0]'; 
    n(:,4)=[0 0 0]'; 
    
    
    r2=[cos(xbase(i-1)) -sin(xbase(i-1)) 0;sin(xbase(i-1)) cos(xbase(i-1)) 0; 0 0 1]; %set R to i+1
    fnew=r2*fin(:,i)+jointforce(:,i);
    fin(:,i-1)=[fnew];
    c1=r2*n(:,i); 
    c2=N(:,i);
    c3a=r2*fin(:,i);
    c3b=cross(pbase(:,i-1),c3a);
    c4=cross(sbase(:,i-1),jointforce(:,i));
    ninw=c1+c2+c3b+c4;
    n(:,i-1)=[ninw];
    torque=[n(3,[1:3])];
    
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mass, distancecom, masscom] = cylcoms (r1,r2,len,density, exoforce)
%************** MATLAB "cylcoms function (Olivia Ball) *************
% Calculate the mass of a cylinder with varying bases to create a new COM
% when given outside exoskeletal forces 

% location of the center of mass 
% SYNTAX      [mass, distancecom, masscom] = cylcoms (r1,r2,len,density,exoforce)   
% INPUTS:       r1          Radius of one of the end point 
%               r2          Radius of the other endpoint 
%               length      Length of the segment of interest
%               density     Desity of the object of interest
%               exoforce   Forces on the exoskeleton
%             
% OUTPUT:      mass         The original mass of the cylinder 
%              distancecom	Distance from the proximal end to the COM
%              masscom      The mass at the COM
% CALLS:       
% CALLED BY:	
%~~~~~~~~~~~~~~~~~~~~~~ Begin Program: ~~~~~~~~~~~~~~~~~~~~~~~~~~

mass=(((density*pi*len)/3)*(r1^2+r2*r1+r2^2)); % Mass calculation 
zplane=((3*r1^2+2*r1*r2+r2^2)/(r1^2+r1*r2+r2^2))*(len/4); % Zplane
distancecom=((mass*zplane+exoforce*0)/(mass+exoforce)); % Distance of COM from proximal end 
masscom=(exoforce+mass); % Mass for new COM 

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ixx, iyy, izz] = momentinter(r1,r2,distancecom,masscom)  
%************** MATLAB "mementinter function (Olivia Ball) *************
% Calulate the moment of interia for each of the planes

% location of the center of mass 
% SYNTAX      [ixx, iyy, izz] = momentinteria(r1,r2,distancecom,masscom) 
% INPUTS:       r1          Radius of one of the end point 
%               r2          Radius of the other endpoint 
%               distancecom	Distance from the proximal end to the COM
%               masscom      The mass at the COM
%             
% OUTPUT:      ixx         Moment of inertia in the x-direction
%              iyy         Moment of inertia in the y-direction 
%              izz         Moment of inertia in the z-direction
%              
% CALLS:       
% CALLED BY:	
%~~~~~~~~~~~~~~~~~~~~~~ Begin Program: ~~~~~~~~~~~~~~~~~~~~~~~~~~

rnew= ((r1+r2)/2); %radius for cylinder calculation - assumption
ixx=masscom*((distancecom^2/12)+(rnew^2/4)); %kg*m^2, moment of inertia for upper leg in x direction
iyy=masscom*((distancecom^2/12)+(rnew^2/4)); %kg*m^2, moment of inertia for upper leg in y direction
izz= (masscom*rnew^2)/2; %kg*m^2, upper leg

end
