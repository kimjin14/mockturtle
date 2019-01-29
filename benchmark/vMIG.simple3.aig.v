// Test
module top (pi0, pi1, pi2, pi3, pi4, pi5, pi6, pi7, pi8, pi9, pi10, pi11, pi12, pi13, pi14, pi15, pi16, pi17, po0);
input pi0, pi1, pi2, pi3, pi4, pi5, pi6, pi7, pi8, pi9, pi10, pi11, pi12, pi13, pi14, pi15, pi16, pi17;
output po0;
wire w0, w1, w2, w3, w4, w5, w6, w7, w8, w9;
assign w0 = (pi0 & pi1) | (pi0 & pi2) | (pi1 & pi2);
assign w1 = (pi3 & pi4) | (pi3 & pi5) | (pi4 & pi5);
assign w2 = (pi6 & pi7) | (pi6 & pi8) | (pi7 & pi8);
assign w3 = (pi9 & pi10) | (pi9 & pi11) | (pi10 & pi11);
assign w4 = (pi12 & pi13) | (pi12 & pi14) | (pi13 & pi14);
assign w5 = (pi15 & pi16) | (pi15 & pi17) | (pi16 & pi17);
assign w9 = (pi3 & pi4) | (pi3 & pi8) | (pi4 & pi8);
assign w6 = (w0 & w1) | (w0 & w2) | (w1 & w2);
assign w7 = (w6 & w2) | (w6 & w9) | (w2 & w9);
assign w8 = (w7 & w3) | (w7 & w4) | (w3 & w4);
assign po0 = (w8 & w4) | (w8 & w5) | (w4 & w5);
endmodule
