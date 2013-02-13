package com.github.mlc.jgeographic;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class MathTest {
    @Test
    public void testPi() {
        assertEquals(3.14159, java.lang.Math.PI, 1e-5);
        assertEquals(9.86960, Math.sq(java.lang.Math.PI), 1e-5);
    }

    @Test
    public void testAtanh() {
        assertEquals(1.0986123, Math.atanh(0.8), 1e-7);
    }

    @Test
    public void testAsinh() {
        assertEquals(0.0, Math.asinh(0.0), 1e-40);
        assertEquals(0.88137359, Math.asinh(1.0), 1e-8);
        assertEquals(java.lang.Math.PI / 4, Math.asinh(0.86867096), 1e-7);
    }
}
