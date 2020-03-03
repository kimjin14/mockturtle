//Written by the Majority Logic Package Thu Apr 30 13:27:53 2015
module top (pi0, pi1, pi2, pi3, pi4, pi5, pi6, pi7, pi8, pi9, pi10, pi11, pi12, po0);
input pi0, pi1, pi2, pi3, pi4, pi5, pi6, pi7, pi8, pi9, pi10, pi11, pi12;
output po0;
wire w0, w1, w2, w3, w4;

assign w0 = (pi0 & pi1) | (pi0 & pi2) | (pi1 & pi2);
assign w1 = (pi3 & pi4) | (pi3 &  w0) | (pi4 &  w0);
assign w2 = (pi5 & pi6) | (pi5 &  w1) | (pi6 &  w1);
assign w3 = (~pi7 & pi8) | (~pi7 &  ~w2) | (pi8 &  ~w2);
assign w4 = (~pi9 & pi10) | (~pi9 & w3) | (pi10 & w3);
assign po0 = (pi11 & pi12) | (pi11 & w4) | (pi12 & w4);

endmodule
