%%%%% Run desiredtorques.m FIRST and do not clear the workspace %%%%%%%
x_off = (lu*100)+(ll*100);
y_off = 0;

for j = 1:26
    pos_x(j) = x_off-((j-1)/100)^2*80;
    pos_y(j) = y_off+((j-1)/100)^2*80;
    rot(j) = ((j-1)/100)^2*80*pi/180;
    timec(j) = (j-1)/100;
end
x_current = pos_x(26);
y_current = pos_y(26);
rot_current = rot(26);

for j = 27:76
    pos_x(j) = x_current-(j-26)/100*40;
    pos_y(j) = y_current+(j-26)/100*40;
    rot(j) = rot_current+(j-26)/100*40*pi/180;
    timec(j) = (j-1)/100;
end

x_current = pos_x(76);
y_current = pos_y(76);
rot_current = rot(76);

for j = 77:101
    pos_x(j) = x_current+((j-76)/100)^2*80-(j-76)/100*40;
    pos_y(j) = y_current-((j-76)/100)^2*80+(j-76)/100*40;
    rot(j) = rot_current+(-((j-76)/100)^2*80+(j-76)/100*40)*pi/180;
    timec(j) = (j-1)/100;
end

output = [timec',pos_x',pos_y',rot'];

% Create the text file 
fileID = fopen('trajoutput.txt','w');
fprintf(fileID,'%5f %5f %5f %5f\n',[output].');
fclose(fileID);
