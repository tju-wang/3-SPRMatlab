

function Fu =solve_fun(x)
    global gX0 gY0 gZ0;
    s = 62;  %测量得 58.5mm + 6.5/2
    a = 41.56; %动平台外接圆半径
    b = 80;
    Fu(1) =( gX0 - (2*(gZ0)*(cos(x(1))*sin(x(2))*(cos(x(1))*cos(x(3)) - cos(x(2))*sin(x(1))*sin(x(3))) + sin(x(1))*sin(x(2))*(cos(x(3))*sin(x(1)) + cos(x(1))*cos(x(2))*sin(x(3)))) - b*sin(x(1))*sin(x(2))*(cos(x(2)) - 3*cos(x(1))*cos(x(3)) + 3*cos(x(2))*sin(x(1))*sin(x(3))))/(2*cos(x(2))*(cos(x(1))*cos(x(3)) - cos(x(2))*sin(x(1))*sin(x(3))) - 2*sin(x(1))*sin(x(2))^2*sin(x(3))));
    Fu(2) =( gY0 - -(2*(gZ0)*(cos(x(2))*(cos(x(3))*sin(x(1)) + cos(x(1))*cos(x(2))*sin(x(3))) + cos(x(1))*sin(x(2))^2*sin(x(3))) - b*cos(x(2))*(cos(x(2)) - cos(x(1))*cos(x(3)) + cos(x(2))*sin(x(1))*sin(x(3))) + 2*b*sin(x(1))*sin(x(2))^2*sin(x(3)))/(2*cos(x(2))*(cos(x(1))*cos(x(3)) - cos(x(2))*sin(x(1))*sin(x(3))) - 2*sin(x(1))*sin(x(2))^2*sin(x(3))));
    Fu(3) =(x(1) - x(3));
end