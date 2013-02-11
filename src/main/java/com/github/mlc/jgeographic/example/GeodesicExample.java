/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic.example;

import com.github.mlc.jgeographic.Constants;
import com.github.mlc.jgeographic.Geodesic;
import com.github.mlc.jgeographic.InversePosition;
import com.github.mlc.jgeographic.LatLon;

public class GeodesicExample {
    public static void main(String[] args) {
        Geodesic geod = new Geodesic(Constants.WGS84_a, Constants.WGS84_f);
        // alternately: Geodesic geod = Geodesic.WGS84;
        {
            // Sample direct calculation, travelling about NE from JFK
            double lat1 = 40.6, lon1 = -73.8, s12 = 5.5e6, azi1 = 51;
            LatLon position = new LatLon();
            geod.direct(lat1, lon1, azi1, s12, position);
            System.out.println(position.lat + " " + position.lon);
        }
        {
            // Sample inverse calculation, JFK to LHR
            double
                    lat1 = 40.6, lon1 = -73.8, // JFK Airport
                    lat2 = 51.6, lon2 = -0.5;  // LHR Airport
            InversePosition inversePosition = new InversePosition();
            geod.genInverse(lat1, lon1, lat2, lon2, Geodesic.Mask.DISTANCE, inversePosition);
            System.out.println(inversePosition.s12);
        }
    }
}
