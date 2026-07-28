// Knuth's model of a fair die using only fair coins
dtmc

module coin
    side : [0..1] init 0;

    [toss] true -> 0.5 : (side'=0) + 0.5 : (side'=1);
endmodule

module die
	// local state
	s : [0..7] init 0;
	// value of the dice
	d : [0..6] init 0;
    	// phase
    	flip : [0..2] init 1;
	
    [toss] flip=1 -> 1: (flip'=0);

	[move] s=0 & side=0 & flip=0 -> 1 : (s'=1) & (flip'=1);
	[move] s=1 & side=0 & flip=0 -> 1 : (s'=3) & (flip'=1);
	[move] s=2 & side=0 & flip=0 -> 1 : (s'=5) & (flip'=1);
	[move] s=3 & side=0 & flip=0 -> 1 : (s'=1) & (flip'=1);
	[move] s=4 & side=0 & flip=0 -> 1 : (s'=7) & (d'=2) & (flip'=2);
	[move] s=5 & side=0 & flip=0 -> 1 : (s'=7) & (d'=4) & (flip'=2);
	[move] s=6 & side=0 & flip=0 -> 1 : (s'=2) & (flip'=1);
	
    [move] s=7 -> (s'=7);

    [move] s=0 & side=1 & flip=0 -> 1 : (s'=2) & (flip'=1);
	[move] s=1 & side=1 & flip=0 -> 1 : (s'=4) & (flip'=1);
	[move] s=2 & side=1 & flip=0 -> 1 : (s'=6) & (flip'=1);
	[move] s=3 & side=1 & flip=0 -> 1 : (s'=7) & (d'=1) & (flip'=2);
	[move] s=4 & side=1 & flip=0 -> 1 : (s'=7) & (d'=3) & (flip'=2);
	[move] s=5 & side=1 & flip=0 -> 1 : (s'=7) & (d'=5) & (flip'=2);
	[move] s=6 & side=1 & flip=0 -> 1 : (s'=7) & (d'=6) & (flip'=2);
	
endmodule

rewards "coin_flips"
	[toss] s<7 : 1;
endrewards


label "target" = s=7 & d=3;
label "flip_heads" = flip=0 & side=0;
label "flip_tails" = flip=0 & side=1;