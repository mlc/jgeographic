package com.github.mlc.jgeographic;

import org.junit.Test;
import static org.junit.Assert.*;

public class DMSTest {
    @Test
    public void canDecodeLatitude() {
        String dms = "30d14'45.6\"S";
        DecodeResult dr = DMS.decode(dms);
        assertEquals(DMS.Flag.LATITUDE, dr.flag);
        assertEquals(-30.246, dr.n, 1e-10);
    }

    @Test
    public void canDecodeLongitude() {
        String dms = "74°0′23.04″W";
        DecodeResult dr = DMS.decode(dms);
        assertEquals(DMS.Flag.LONGITUDE, dr.flag);
        assertEquals(-74.0064, dr.n, 1e-10);
    }

    @Test
    public void canEncode() {
        double angle = -30.245715;
        assertEquals("30d14'44.57\"S", DMS.encode(angle, 6, DMS.Flag.LATITUDE));
    }
}
