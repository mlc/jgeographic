/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic;

import static java.lang.Math.*;

public class Math {
    private Math() { }

    public static final double degree = PI / 180.0;
    public static double sq(double x) { return x * x; }
    public static int sq(int x) { return x * x; }

    public static double asinh(double x) {
        double y = abs(x);
        y = log1p(y * (1 + y /(hypot(1.0, y) + 1)));
        return x < 0 ? -y : y;
    }

    public static double atanh(double x) {
        double y = abs(x);
        y = log1p(2 * y/(1 - y))/2;
        return x < 0 ? -y : y;
    }

    public static DoublePair sum(double u, double v) {
        double s = u + v;
        double up = s - v;
        double vpp = s - up;
        up -= u;
        vpp -= v;
        double t = -(up + vpp);
        return new DoublePair(s, t);
    }

    /**
     * Normalize an angle (restricted input range).
     *
     * @param x the angle in degrees, must be in [&minus;540&deg;, 540&deg;)
     * @return the angle reduced to the range [&minus;180&deg;, 180&deg;]
     */
    public static double angNormalize(double x) {
        return x >= 180 ? x - 360 : (x < -180 ? x + 360 : x);
    }

    /**
     * Normalize an angle (arbitrary input range).
     *
     * @param x the angle in degrees (unrestricted range)
     * @return the angle reduced to the range [&minus;180&deg;, 180&deg;]
     */
    public static double angNormalize2(double x) {
        return angNormalize(x % 360.0);
    }

    /**
     * Difference of two angles reduced to [&minus;180&deg;, 180&deg;]
     *
     * <i>x</i> and <i>y</i> must both lie in [&minus;180&deg;, 180&deg;].  The result
     * is equivalent to computing the difference exactly, reducing it to
     * (&minus;180&deg;, 180&deg;] and rounding the result.  Note that this
     * prescription allows &minus;180&deg; to be returned (e.g., if <i>x</i> is
     * tiny and negative and <i>y</i> = 180&deg;).

     * @param x the first angle, in degrees
     * @param y the second angle, in degrees
     * @return <i>y</i> &minus; <i>x</i>, reduced to the range [&minus;180&deg;,
     *   180&deg;].
     */
    public static double angDiff(double x, double y) {
        // inline sum() here so we don't have to create an object for its return
        double u = -x;
        double d = u + y;
        double up = d - y;
        double vpp = d - up;
        up -= u;
        vpp -= y;
        double t = -(up + vpp);
        if ((d - 180.0) + t > 0.0) // y - x > 180
            d -= 360.0;            // exact
        else if ((d + 180.0) + t <= 0.0) // y - x <= -180
            d += 360.0;            // exact
        return d + t;
    }
    public static boolean isFinite(double d) {
        return !Double.isInfinite(d);
    }
}
