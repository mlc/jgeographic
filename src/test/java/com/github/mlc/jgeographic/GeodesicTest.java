package com.github.mlc.jgeographic;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class GeodesicTest {
    private Geodesic geod;

    @Before
    public void setGeod() {
        geod = new Geodesic(Constants.WGS84_a, Constants.WGS84_f);
    }

    @Test
    public void testDirect() {
        double lat1 = 40.6, lon1 = -73.8, s12 = 5.5e6, azi1 = 51;
        LatLon position = new LatLon();
        geod.direct(lat1, lon1, azi1, s12, position);
        assertEquals(51.8846, position.lat, 1e-4);
        assertEquals(-1.14117, position.lon, 1e-4);
    }

    @Test
    public void testInverse() {
        double
                lat1 = 40.6, lon1 = -73.8, // JFK Airport
                lat2 = 51.6, lon2 = -0.5;  // LHR Airport
        InversePosition inversePosition = new InversePosition();
        geod.genInverse(lat1, lon1, lat2, lon2, Geodesic.Mask.DISTANCE, inversePosition);
        assertEquals(5.55176e+06, inversePosition.s12, 10);
    }
}
