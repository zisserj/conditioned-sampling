// Based on GRID WORLD MODEL OF ROBOT AND JANITOR
// Hakan Younes/gxn/dxp 04/05/04
dtmc

// CONSTANTS
const int xn = 6; // size of the grid
const int yn = 5;

// the following formulae return 1 or 0 depending on whether
// the robot can move in that direction or not
formula right = min(1,max(xn-x1,0));
formula left = min(1,max(x1-1,0));
formula up = min(1,max(yn-y1,0));
formula down = min(1,max(y1-1,0));

// total number of moves the robot randomly picks from
formula moves = right+left+up+down;
module robot
	x1 : [1..xn] init 1; // x position of robot (bottom left)
	y1 : [1..yn] init 1; // y position of robot
	
	[move] true    -> up/moves:(y1'=y1+1) + down/moves:(y1'=y1-1) + left/moves:(x1'=x1-1) + right/moves:(x1'=x1+1); 
	[bump_down] (on_db1|on_db2|on_db3|on_db4) & (down=1) -> 1: (y1'=y1-1);
    [bump_up] (on_ub5|on_ub6|on_ub7|on_ub8) & (up=1) -> 1: (y1'=y1+1);
endmodule
module down_bumper1
	cx1 : [1..xn] init 2;
	cy1 : [1..yn] init 2;
	[bump_down] true -> 1: true;
	[move] !on_db1 -> 1: true; 
endmodule
formula on_db1 = cx1 = x1 & cy1 = y1;

module down_bumper2
	cx2 : [1..xn] init 3;
	cy2 : [1..yn] init 2;
	[bump_down] true -> 1: true;
	[move] !on_db2 -> 1: true; 
endmodule
formula on_db2 = cx2 = x1 & cy2 = y1;

module down_bumper3
	cx3 : [1..xn] init 4;
	cy3 : [1..yn] init 2;
	[bump_down] true -> 1: true;
	[move] !on_db3 -> 1: true; 
endmodule
formula on_db3 = cx3 = x1 & cy3 = y1;

module down_bumper4
	cx4 : [1..xn] init 5;
	cy4 : [1..yn] init 2;
	[bump_down] true -> 1: true;
	[move] !on_db4 -> 1: true; 
endmodule
formula on_db4 = cx4 = x1 & cy4 = y1;
module up_bumper5
	cx5 : [1..xn] init 2;
	cy5 : [1..yn] init 4;
	[bump_up] true -> 1: true;
	[move] !on_ub5 -> 1: true; 
endmodule
formula on_ub5 = cx5 = x1 & cy5 = y1;

module up_bumper6
	cx6 : [1..xn] init 3;
	cy6 : [1..yn] init 4;
	[bump_up] true -> 1: true;
	[move] !on_ub6 -> 1: true; 
endmodule
formula on_ub6 = cx6 = x1 & cy6 = y1;

module up_bumper7
	cx7 : [1..xn] init 4;
	cy7 : [1..yn] init 4;
	[bump_up] true -> 1: true;
	[move] !on_ub7 -> 1: true; 
endmodule
formula on_ub7 = cx7 = x1 & cy7 = y1;

module up_bumper8
	cx8 : [1..xn] init 5;
	cy8 : [1..yn] init 4;
	[bump_up] true -> 1: true;
	[move] !on_ub8 -> 1: true; 
endmodule
formula on_ub8 = cx8 = x1 & cy8 = y1;
label "target" = x1=6 & y1=5;