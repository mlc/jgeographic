/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic;

public class Constants {
    private Constants() { }
    /** The number of radians in a degree. */
    public static final double degree = Math.degree;
    /** The number of radians in an arcminute. */
    public static final double arcminute = degree / 60.0;
    /** The number of radians in an arcsecond. */
    public static final double arcsecond = degree / 3600.0;

    public static final double meter = 1.0,
            kilometer = 1000.0 * meter,
            nauticalMile = 1852 * meter,
            squareMeter = meter * meter,
            hectatre = 10000 * squareMeter,
            squareKilometer = kilometer * kilometer,
            squareNauticalMile = nauticalMile * nauticalMile,
            foot = 0.0254 * 12 * meter,
            yard = 3 * foot,
            fathom = 2 * yard,
            chain = 22 * yard,
            furlong = 10 * chain,
            mile = 8 * furlong,
            acre = chain * furlong,
            squareMile = mile * mile,
            surveyFoot = 1200.0 / 3937.0 * meter;

    public static final double WGS84_a = 6378137 * meter;
    public static final double WGS84_f = 1/298.257223563;
    public static final double WGS84_GM = 3.986004418E14;
    public static final double WGS84_omega = 729211.5;
    public static final double WGS84_r = 1/WGS84_f;

    public static final double GRS80_a = 6378137 * meter;
    public static final double GRS80_GM = 3.986005E14;
    public static final double GRS80_omega = 729211.5;
    public static final double GRS80_J2 = 0.00108263;

    public static final double UTM_k0 = 0.9996;
    public static final double UPS_k0 = 0.994;
}
