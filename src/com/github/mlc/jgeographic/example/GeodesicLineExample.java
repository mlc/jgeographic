/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic.example;

import java.text.NumberFormat;
import java.util.Locale;
import com.github.mlc.jgeographic.Geodesic;
import com.github.mlc.jgeographic.GeodesicLine;
import com.github.mlc.jgeographic.InversePosition;
import com.github.mlc.jgeographic.LatLon;

public class GeodesicLineExample {
    public static void main(String[] args) {
        // Print waypoints between JFK and SIN
        Geodesic geod = Geodesic.WGS84;
        double
                lat1 = 40.640, lon1 = -73.779, // JFK
                lat2 =  1.359, lon2 = 103.989; // SIN
        double s12, azi1, azi2;
        InversePosition inverse = new InversePosition();
        double a12 = geod.genInverse(lat1, lon1, lat2, lon2, Geodesic.Mask.AZIMUTH | Geodesic.Mask.DISTANCE, inverse);
        s12 = inverse.s12;
        azi1 = inverse.azi1;
        azi2 = inverse.azi2;
        GeodesicLine line = geod.line(lat1, lon1, azi1);
        double ds = 500e3;
        int num = (int)Math.ceil(s12 / ds);
        NumberFormat nf = NumberFormat.getInstance(Locale.US);
        nf.setMinimumFractionDigits(3);
        nf.setMaximumFractionDigits(3);
        {
            // use intervals of equal length
            double ds_s = s12 / num;
            for (int i = 0; i <= num; ++i) {
                LatLon position = line.position(i * ds_s);
                System.out.println(i + " " + nf.format(position.lat) + " " + nf.format(position.lon));
            }
        }
        {
            // slightly faster: use intervals of equal arclength
            double da = a12 / num;
            for (int i = 0; i <= num; ++i) {
                LatLon position = line.arcPosition(i * da);
                System.out.println(i + " " + nf.format(position.lat) + " " + nf.format(position.lon));
            }
        }
    }
}
