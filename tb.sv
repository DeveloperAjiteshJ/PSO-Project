`timescale 1ns/1ps

module tb_pso_rc_tuner;

  // Override parameters for this test
  pso_rc_tuner #(
    .N_PARTICLES(15),
    .N_ITER     (50),
    .F_DESIRED  (12345.0),
    .R_MIN      (1.0e2),
    .R_MAX      (1.0e8),
    .C_MIN      (100e-12),
    .C_MAX      (1.0e-6),
    .C1         (1.6),
    .C2         (1.6),
    .W_START    (0.7),
    .W_END      (0.7)
  ) dut ();

  // Safety timeout: 100 ms at 1 ns timescale
  initial begin
    #100_000_000;
    $display("Timeout: stopping simulation");
    $finish;
  end

endmodule
