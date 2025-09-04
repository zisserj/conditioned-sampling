
# Remember model is 1-idx

gridsize = (6, 5)
down_bumpers = [(i, 2) for i in range(2, 6)]
up_bumpers = [(i, 4) for i in range(2, 6)]

robot_init = (1,1)
robot_goal = (6, 5)

filename = 'robot_auto.pm'

program_header = '''// Based on GRID WORLD MODEL OF ROBOT AND JANITOR
// Hakan Younes/gxn/dxp 04/05/04
dtmc

// CONSTANTS
const int xn = {}; // size of the grid
const int yn = {};

// the following formulae return 1 or 0 depending on whether
// the robot can move in that direction or not
formula right = min(1,max(xn-x1,0));
formula left = min(1,max(x1-1,0));
formula up = min(1,max(yn-y1,0));
formula down = min(1,max(y1-1,0));

// total number of moves the robot randomly picks from
formula moves = right+left+up+down;
'''.format(*gridsize)

program_robot = '''module robot
	x1 : [1..xn] init 1; // x position of robot (bottom left)
	y1 : [1..yn] init 1; // y position of robot
	
	[move] true    -> up/moves:(y1'=y1+1) + down/moves:(y1'=y1-1) + left/moves:(x1'=x1-1) + right/moves:(x1'=x1+1); 
	[bump_down] ({}) & (down=1) -> 1: (y1'=y1-1);
    [bump_up] ({}) & (up=1) -> 1: (y1'=y1+1);
endmodule
'''.format('|'.join([f'on_db{i+1}' for i in range(len(down_bumpers))]),
           '|'.join([f'on_ub{i+1+len(down_bumpers)}' for i in range(len(up_bumpers))]))

db_module = lambda i, pos: f'''module down_bumper{i}
	cx{i} : [1..xn] init {pos[0]};
	cy{i} : [1..yn] init {pos[1]};
	[bump_down] true -> 1: true;
	[move] !on_db{i} -> 1: true; 
endmodule
formula on_db{i} = cx{i} = x1 & cy{i} = y1;
'''

ub_module = lambda i, pos: f'''module up_bumper{i}
	cx{i} : [1..xn] init {pos[0]};
	cy{i} : [1..yn] init {pos[1]};
	[bump_up] true -> 1: true;
	[move] !on_ub{i} -> 1: true; 
endmodule
formula on_ub{i} = cx{i} = x1 & cy{i} = y1;
'''

target = 'label "target" = x1={} & y1={};'.format(*robot_goal)

program_db = '\n'.join([db_module(i, pos) for i, pos in enumerate(down_bumpers, 1)])
program_ub = '\n'.join([ub_module(i, pos) for i, pos in enumerate(up_bumpers, 1+len(down_bumpers))])

with open(filename, 'w+') as f:
    f.writelines([program_header, program_robot, program_db, program_ub, target])
