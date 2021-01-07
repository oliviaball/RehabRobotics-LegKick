%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond = textread('trajoutput.txt'); % Read values created for given leg
samp_interval=0.01; %Define the interval period 
time=[samp_interval*(0:100)]'; %Create the correct time outputs
x = [ cond(:,2)' cond(:,3)' cond(:,4)'];
deriv = deriv3pt(x, samp_interval);
x_deriv=deriv(1:101);
y_deriv=deriv(102:202);
theta=deriv(203:303);

%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------ Will not automatically update ------------------
format short

syms theta1 theta2 theta3

% Define the f1 and f2 so can take partial derivative for Jacobian
f1=lu*cos(theta1)+ll*((cos(theta1)*cos(theta2)-(sin(theta1)*sin(theta2))));
f2=lu*sin(theta1)+ll*((cos(theta1)*sin(theta2))+cos(theta2)*sin(theta1));
f3=0;

% Take the partial derivative
f11=diff(f1,theta1);
f12=diff(f1,theta2);
f13=diff(f1,theta3);
f21=diff(f2,theta1);
f22=diff(f2,theta2);
f23=diff(f2,theta3);
%----------------------------------------------------------------------

% Create a locker for each value to be held in
thetadot_final=[];
thetadot1=[];
thetadot2=[];
thetadot3=[];
x1 = zeros([1,length(x_deriv)]);
x2 = zeros([1,length(y_deriv)]);
x3 = zeros([1,length(theta)]);

for j=2:101
    
    %Set the initial conditions
    x1(1)=0.05; % Theta 1 (radians)
    x2(1)=-0.05; % Theta 2 ( radians)
    x3(1)=0; % Theta 3 (radians)
    
    % No initial velocity assumed %
    thetadot1(:,1)=0;
    thetadot2(:,1)=0;
    thetadot3(:,1)=0;
    
    % Pull the derivatives from derivpt3 to multiply by the jacobian
    a2=[x_deriv(j-1); y_deriv(j-1); theta(j-1)]; %3X1 matrix to get thetadot values
   
    %define x1 and x2 instead of theta1 and theta2
    top=[-36*sin(x1(:,j-1))-52*cos(x1(:,j-1))*sin(x2(:,j-1))-52*cos(x2(:,j-1))*sin(x1(:,j-1)) -52*cos(x1(:,j-1))*sin(x2(:,j-1))-52*cos(x2(:,j-1))*sin(x1(:,j-1)) 0];
    middle=[36*cos(x1(:,j-1))+52*cos(x1(:,j-1))*cos(x2(:,j-1))-52*sin(x1(:,j-1))*sin(x2(:,j-1)) 52*cos(x1(:,j-1))*cos(x2(:,j-1))-52*sin(x1(:,j-1))*sin(x2(:,j-1)) 0];
    bottom=[1 1 1];
    J0=[top;middle;bottom]; %define the jacobian in terms of x1 and x2 so can change the theta values each rough
    thetadot=inv(J0)*a2; %inverse matrix * derivatives to get thetadot
    thetadot_final(:,j)=[thetadot]; %create a matrix for each of the thetadots at each time stamp
    
    %pull just thetadot1, thetadot2, thetadot3
    thetadot1(:,j)=[thetadot(1,:)];
    thetadot2(:,j)=[thetadot(2,:)];
    thetadot3(:,j)=[thetadot(3,:)]; 
    
    %take the interval at each time interval to get theta1, theta2, theta3
    %1/2*delta(time) where delta(time)=0.01
    theta1=(samp_interval*0.5)*(thetadot1(:,j-1)+thetadot1(:,j)); %thetadot1(j-1)+thetadot1(j) is two lengths then "height" is 0.01 and times 1/2
    x1(:,j)=x1(:,j-1)+theta1; %old theta value+new theta value
   
    theta2=((samp_interval*0.5)*(thetadot2(:,j-1)+thetadot2(:,j))); %thetadot2(j-1)+thetadot2(j) is two lengths then "height" is 0.01 and times 1/2
    x2(:,j)=x2(:,j-1)+theta2; %old theta value+new theta value
    
    theta3=(samp_interval*0.5)*(thetadot3(:,j-1)+thetadot3(:,j)); %thetadot3(j-1)+thetadot3(j) is two lengths then "height" is 0.01 and times 1/2
    x3(:,j)=x3(:,j-1)+theta3; %old theta value+new theta value
    
   sumtheta=[x1+x2+x3]; %sum of all the theta values

end

finaljoints = [time,x1.',x2.', x3.']

filejoints = fopen('lastname_jointangles.txt','w');
fprintf(filejoints,'%5f %5f %5f %5f\n',[finaljoints].');
fclose(filejoints);
%%
%%%%%% All the plots to make sure that everything looks right%%%%%%
figure (1)
subplot(3,1,1)
plot(time,cond(:,2)')
title('Desired Moment for the Desired Trajectory')
xlabel('Time (Seconds)')
ylabel('Magnitude of X (cm)')
subplot(3,1,2)
plot(time,cond(:,3)')
title('Y Movement for the Desired Trajectory')
xlabel('Time (Seconds)')
ylabel('Magnitude of Y (cm)')
subplot(3,1,3)
plot(time,cond(:,4)')
title('Rotation Values for the Desired Trajectory')
xlabel('Time (Seconds)')
ylabel('Magnitude of Rotation (rad)')

%%%%%%%%% Plot of the Derivatives %%%%%%%%%%%
figure(2)
subplot(3,1,1)

plot(time,x_deriv)
title('Plot of the Derivative of X vs. Time')
xlabel('Time (Seconds)')
ylabel('Velocity in the X-Direction')

subplot(3,1,2)
plot(time,y_deriv)
title('Plot of the Derivative of Y vs. Time')
xlabel('Time (Seconds)')
ylabel('Velocity in the Y-Direction')

subplot(3,1,3)

plot(time,theta)
title('Plot of Derivative of Rotation vs. Time')
xlabel('Time (Seconds)')
ylabel('Velocity of Rotation')

%%%%%%% Plot of Theta Values at each time point %%%%%%%%
figure(3)
subplot(4,1,1)
plot(time,x1)
title('Plot of Theta 1 vs. Time')
xlabel('Time (Seconds)')
ylabel('Angle Theta 1')

subplot(4,1,2)
plot(time,x2)
title('Plot of Theta 2 vs. Time')
xlabel('Time (Seconds)')
ylabel('Angle Theta 2')

subplot(4,1,3)
plot(time,x3)
title('Plot of Theta 3 vs. Time')
xlabel('Time (Seconds)')
ylabel('Angle Theta 3')

subplot(4,1,4)
plot(time,sumtheta)
title('Sum of Theta Values')
xlabel('Time (Seconds)')
ylabel('Angle Sum of Theta')


