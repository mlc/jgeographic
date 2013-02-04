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
import static com.github.mlc.jgeographic.Math.*;

public class Geodesic {
    // Underflow guard
    static final double tiny = sqrt(Double.MIN_VALUE);
    static final double tol0 = 2.22044604925031308085e-16;
    static final double tol1 = 200 * tol0;
    static final double tol2 = sqrt(tol0);

    // Check on bisection interval
    static final double tolb = tol0 * tol2;
    static final double xthresh = 1000 * tol2;

    static final int GEODESIC_ORDER = 6;
    static final int nA1_ =  GEODESIC_ORDER;
    static final int nC1_ =  GEODESIC_ORDER;
    static final int nC1p_ = GEODESIC_ORDER;
    static final int nA2_ =  GEODESIC_ORDER;
    static final int nC2_ =  GEODESIC_ORDER;
    static final int nA3_ =  GEODESIC_ORDER;
    static final int nA3x_ = nA3_;
    static final int nC3_ =  GEODESIC_ORDER;
    static final int nC3x_ = (nC3_ * (nC3_ - 1)) / 2;
    static final int nC4_ =  GEODESIC_ORDER;
    static final int nC4x_ = (nC4_ * (nC4_ + 1)) / 2;
    static final int maxit1_ = 20;
    static final int maxit2_ = 93;

    static final int CAP_NONE = 0,
            CAP_C1   = 1<<0,
            CAP_C1p  = 1<<1,
            CAP_C2   = 1<<2,
            CAP_C3   = 1<<3,
            CAP_C4   = 1<<4,
            CAP_ALL  = 0x1F,
            OUT_ALL  = 0x7F80;

    static double angRound(double x) {
        final double z = 1/16.0;
        double y = abs(x);
        y = y < z ? z - (z - y) : y;
        return x < 0 ? -y : y;
    }

    static void sinCosNorm(DoublePair sinxcosx) {
        double r = hypot(sinxcosx.s, sinxcosx.c);
        sinxcosx.s /= r;
        sinxcosx.c /= r;
    }

    final double _a, _f, _f1, _e2, _ep2, _n, _b, _c2, _etol2;
    final double[] _A3x = new double[nA3x_], _C3x = new double[nC3x_], _C4x = new double[nC4x_];

    public static class Mask {
        private Mask() { }
        public static final int NONE = 0;
        public static final int LATITUDE =      1<<7  | CAP_NONE;
        public static final int LONGITUDE =     1<<8  | CAP_C3;
        public static final int AZIMUTH =       1<<9  | CAP_NONE;
        public static final int DISTANCE =      1<<10 | CAP_C1;
        public static final int DISTANCE_IN =   1<<11 | CAP_C1 | CAP_C1p;
        public static final int REDUCEDLENGTH = 1<<12 | CAP_C1 | CAP_C2;
        public static final int GEODESICSCALE = 1<<13 | CAP_C1 | CAP_C2;
        public static final int AREA          = 1<<14 | CAP_C4;
        public static final int ALL           = OUT_ALL | CAP_ALL;
    }

    public Geodesic(double a, double f) {
        _a = a;
        _f = f <= 1 ? f : 1/f;
        _f1 = 1- _f;
        _e2 = _f * (2 - _f);
        _ep2 = _e2 / sq(_f1);
        _n = _f / (2 - _f);
        _b = _a * _f1;
        _c2 = (sq(_a) + sq(_b) *
           (_e2 == 0 ? 1 :
            (_e2 > 0 ? atanh(sqrt(_e2)) : atan(sqrt(-_e2))) /
            sqrt(abs(_e2))))/2; // authalic radius squared
        _etol2 = 0.01 * tol2 / max(0.1, sqrt(abs(_e2)));
        if (!(isFinite(_a) && _a > 0))
            throw new GeographicException("Major radius is not positive");
        if (!(isFinite(_b) && _b > 0))
            throw new GeographicException("Minor radius is not positive");
        A3coeff();
        C3coeff();
        C4coeff();
    }

    public static final Geodesic WGS84 = new Geodesic(Constants.WGS84_a, Constants.WGS84_f);

    public static double sinCosSeries(boolean sinp, double sinx, double cosx, double[] c, int n) {
        // Evaluate
        // y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
        //            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
        // using Clenshaw summation.  N.B. c[0] is unused for sin series
        // Approx operation count = (n + 5) mult and (2 * n + 2) add
        int cptr = (n + (sinp ? 1 : 0));            // Point to one beyond last element
        double ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
                y0 = (n & 1) != 0 ? c[--cptr] : 0, y1 = 0;          // accumulators for sum
        // Now n is even
        n /= 2;
        while (n-- != 0) {
            // Unroll loop x 2, so accumulators return to their original role
            y1 = ar * y0 - y1 + c[--cptr];
            y0 = ar * y1 - y0 + c[--cptr];
        }
        return sinp
                ? 2 * sinx * cosx * y0    // sin(2 * x) * y0
                : cosx * (y0 - y1);       // cos(x) * (y0 - y1)
    }

    public GeodesicLine line(double lat1, double lon1, double azi1, int caps) {
        return new GeodesicLine(this, lat1, lon1, azi1, caps);
    }

    public GeodesicLine line(double lat1, double lon1, double azi1) {
        return new GeodesicLine(this, lat1, lon1, azi1, Mask.ALL);
    }

    public double genDirect(double lat1, double lon1, double azi1,
                            boolean arcmode, double s12_a12, int outmask,
                            Position outPosition) {
        return line(lat1, lon1, azi1, outmask | (arcmode ? Mask.NONE : Mask.DISTANCE_IN)).genPosition(arcmode, s12_a12, outmask, outPosition);
    }
    public double direct(double lat1, double lon1, double azi1, double s12, int outmask, Position outPosition) {
        return genDirect(lat1, lon1, azi1, false, s12, outmask, outPosition);
    }
    public double direct(double lat1, double lon1, double azi1, double s12, LatLon outLatLon) {
        Position outPosition = new Position();
        double arcLength = direct(lat1, lon1, azi1, s12, Mask.LATITUDE | Mask.LONGITUDE, outPosition);
        outLatLon.lat = outPosition.lat2;
        outLatLon.lon = outPosition.lon2;
        return arcLength;
    }

    public double genInverse(double lat1, double lon1, double lat2, double lon2,
                             int outmask,
                             InversePosition out) {
        outmask &= OUT_ALL;
        // Compute longitude difference (AngDiff does this carefully).  Result is
        // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
        // east-going and meridional geodesics.
        double lon12 = angDiff(angNormalize(lon1),
                               angNormalize(lon2));
        // If very close to being on the same half-meridian, then make it so.
        lon12 = angRound(lon12);
        // Make longitude difference positive.
        int lonsign = lon12 >= 0 ? 1 : -1;
        lon12 *= lonsign;
        // If really close to the equator, treat as on equator.
        lat1 = angRound(lat1);
        lat2 = angRound(lat2);
        // Swap points so that point with higher (abs) latitude is point 1
        int swapp = abs(lat1) >= abs(lat2) ? 1 : -1;
        if (swapp < 0) {
            lonsign *= -1;
            double tmp = lat1;
            lat1 = lat2;
            lat2 = tmp;
        }
        // Make lat1 <= 0
        int latsign = lat1 < 0 ? 1 : -1;
        lat1 *= latsign;
        lat2 *= latsign;
        // Now we have
        //
        //     0 <= lon12 <= 180
        //     -90 <= lat1 <= 0
        //     lat1 <= lat2 <= -lat1
        //
        // longsign, swapp, latsign register the transformation to bring the
        // coordinates to this canonical form.  In all cases, 1 means no change was
        // made.  We make these transformations so that there are few cases to
        // check, e.g., on verifying quadrants in atan2.  In addition, this
        // enforces some symmetries in the results returned.

        double phi, sbet1, cbet1, sbet2, cbet2, s12x = Double.NaN, m12x = Double.NaN;

        phi = lat1 * degree;
        // Ensure cbet1 = +epsilon at poles
        DoublePair bet1 = new DoublePair(_f1 * sin(phi), lat1 == -90 ? tiny : cos(phi));
        sinCosNorm(bet1);
        sbet1 = bet1.s;
        cbet1 = bet1.c;

        phi = lat2 * degree;
        // Ensure cbet2 = +epsilon at poles
        DoublePair bet2 = new DoublePair(_f1 * sin(phi), abs(lat2) == 90 ? tiny : cos(phi));
        sinCosNorm(bet2);
        sbet2 = bet2.s;
        cbet2 = bet2.c;

        // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
        // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
        // a better measure.  This logic is used in assigning calp2 in Lambda12.
        // Sometimes these quantities vanish and in that case we force bet2 = +/-
        // bet1 exactly.  An example where is is necessary is the inverse problem
        // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
        // which failed with Visual Studio 10 (Release and Debug)
        if (cbet1 < -sbet1) {
            if (cbet2 == cbet1)
                sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
        } else {
            if (abs(sbet2) == -sbet1)
                cbet2 = cbet1;
        }

        double dn1 = sqrt(1 + _ep2 * sq(sbet1)),
                dn2 = sqrt(1 + _ep2 * sq(sbet2));

        double lam12 = lon12 * degree,
                slam12 = abs(lon12) == 180 ? 0 : sin(lam12),
                clam12 = cos(lam12);      // lon12 == 90 isn't interesting

        double a12 = Double.NaN, sig12, calp1 = Double.NaN, salp1 = Double.NaN, calp2 = Double.NaN, salp2 = Double.NaN;
        // index zero elements of these arrays are unused
        double[] C1a = new double[nC1_ + 1], C2a = new double[nC2_ + 1], C3a = new double[nC3_];
        boolean meridian = lat1 == -90 || slam12 == 0;

        if (meridian) {
            // Endpoints are on a single full meridian, so the geodesic might lie on
            // a meridian.

            calp1 = clam12; salp1 = slam12; // Head to the target longitude
            calp2 = 1; salp2 = 0;           // At the target we're heading north

            double // tan(bet) = tan(sig) * cos(alp)
                    ssig1 = sbet1, csig1 = calp1 * cbet1,
                    ssig2 = sbet2, csig2 = calp2 * cbet2;

            // sig12 = sig2 - sig1
            sig12 = atan2(max(csig1 * ssig2 - ssig1 * csig2, 0.0),
                    csig1 * csig2 + ssig1 * ssig2);
            {
                Lengths l = lengths(_n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                        cbet1, cbet2,
                        (outmask & Mask.GEODESICSCALE) != 0, C1a, C2a);
                s12x = l.s12b;
                m12x = l.m12b;
                out.M12 = l.M12;
                out.M21 = l.M21;
            }
            // Add the check for sig12 since zero length geodesics might yield m12 <
            // 0.  Test case was
            //
            //    echo 20.001 0 20.001 0 | Geod -i
            //
            // In fact, we will have sig12 > pi/2 for meridional geodesic which is
            // not a shortest path.
            if (sig12 < 1 || m12x >= 0) {
                m12x *= _b;
                s12x *= _b;
                a12 = sig12 / degree;
            } else
                // m12 < 0, i.e., prolate and too close to anti-podal
                meridian = false;
        }
        double omg12 = Double.NaN;
        if (!meridian &&
                sbet1 == 0 &&   // and sbet2 == 0
                // Mimic the way Lambda12 works with calp1 = 0
                (_f <= 0 || lam12 <= PI - _f * PI)) {
            // Geodesic runs along equator
            calp1 = calp2 = 0; salp1 = salp2 = 1;
            s12x = _a * lam12;
            sig12 = omg12 = lam12 / _f1;
            m12x = _b * sin(sig12);
            if ((outmask & Mask.GEODESICSCALE) != 0)
                out.M12 = out.M21 = cos(sig12);
            a12 = lon12 / _f1;

        } else if (!meridian) {
            // Now point1 and point2 belong within a hemisphere bounded by a
            // meridian and geodesic is neither meridional or equatorial.

            // Figure a starting point for Newton's method
            {
                DoublePair alp1 = new DoublePair(0.0, 0.0), alp2 = new DoublePair(0.0, 0.0);
                sig12 = inverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                        lam12, alp1, alp2,
                        C1a, C2a);
                salp1 = alp1.s; calp1 = alp1.c;
                salp2 = alp2.s; calp2 = alp2.c;
            }
            if (sig12 >= 0) {
                // Short lines (InverseStart sets salp2, calp2)
                double dnm = (dn1 + dn2) / 2;
                s12x = sig12 * _b * dnm;
                m12x = sq(dnm) * _b * sin(sig12 / dnm);
                if ((outmask & Mask.GEODESICSCALE) != 0)
                    out.M12 = out.M21 = cos(sig12 / dnm);
                a12 = sig12 / degree;
                omg12 = lam12 / (_f1 * dnm);
            } else {
                // Newton's method.  This is a straightforward solution of f(alp1) =
                // lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
                // root in the interval (0, pi) and its derivative is positive at the
                // root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
                // alp1.  During the course of the iteration, a range (alp1a, alp1b) is
                // maintained which brackets the root and with each evaluation of
                // f(alp) the range is shrunk, if possible.  Newton's method is
                // restarted whenever the derivative of f is negative (because the new
                // value of alp1 is then further from the solution) or if the new
                // estimate of alp1 lies outside (0,pi); in this case, the new starting
                // guess is taken to be (alp1a + alp1b) / 2.
                double ssig1 = Double.NaN, csig1 = Double.NaN, ssig2 = Double.NaN, csig2 = Double.NaN, eps = Double.NaN;
                int numit = 0;
                // Bracketing range
                double salp1a = tiny, calp1a = 1, salp1b = tiny, calp1b = -1;
                LambdaResult lr = new LambdaResult();
                for (boolean tripn = false, tripb = false; numit < maxit2_; ++numit) {
                    // the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
                    // WGS84 and random input: mean = 2.85, sd = 0.60
                    double dv;
                    double v = lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                            numit < maxit1_, lr, C1a, C2a, C3a)
                            - lam12;
                    salp2 = lr.salp2;
                    calp2 = lr.calp2;
                    sig12 = lr.sig12;
                    ssig1 = lr.sig1.s;
                    csig1 = lr.sig1.c;
                    ssig2 = lr.sig2.s;
                    csig2 = lr.sig2.c;
                    eps = lr.eps;
                    omg12 = lr.domg12;
                    dv = lr.dlam12;
                    // 2 * tol0 is approximately 1 ulp for a number in [0, pi].
                    // Reversed test to allow escape with NaNs
                    if (tripb || !(abs(v) >= (tripn ? 8 : 2) * tol0)) break;
                    // Update bracketing values
                    if (v > 0 && (numit > maxit1_ || calp1/salp1 > calp1b/salp1b)) {
                        salp1b = salp1; calp1b = calp1;
                    } else if (numit > maxit1_ || calp1/salp1 < calp1a/salp1a) {
                        salp1a = salp1; calp1a = calp1;
                    }
                    if (numit < maxit1_ && dv > 0) {
                        double dalp1 = -v/dv;
                        double sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
                                nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
                        if (nsalp1 > 0 && abs(dalp1) < PI) {
                            DoublePair alp1 = new DoublePair(calp1 * cdalp1 - salp1 * sdalp1, nsalp1);
                            sinCosNorm(alp1);
                            calp1 = alp1.s;
                            salp1 = alp1.c;
                            // In some regimes we don't get quadratic convergence because
                            // slope -> 0.  So use convergence conditions based on epsilon
                            // instead of sqrt(epsilon).
                            tripn = abs(v) <= 16 * tol0;
                            continue;
                        }
                    }
                    // Either dv was not postive or updated value was outside legal
                    // range.  Use the midpoint of the bracket as the next estimate.
                    // This mechanism is not needed for the WGS84 ellipsoid, but it does
                    // catch problems with more eccentric ellipsoids.  Its efficacy is
                    // such for the WGS84 test set with the starting guess set to alp1 =
                    // 90deg:
                    // the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
                    // WGS84 and random input: mean = 4.74, sd = 0.99
                    DoublePair alp1 = new DoublePair((salp1a + salp1b)/2, (calp1a + calp1b)/2);
                    sinCosNorm(alp1);
                    salp1 = alp1.s;
                    calp1 = alp1.c;
                    tripn = false;
                    tripb = (abs(salp1a - salp1) + (calp1a - calp1) < tolb ||
                            abs(salp1 - salp1b) + (calp1 - calp1b) < tolb);
                }
                {
                    Lengths l = lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                            cbet1, cbet2,
                            (outmask & Mask.GEODESICSCALE) != 0, C1a, C2a);
                    s12x = l.s12b;
                    m12x = l.m12b;
                    out.M12 = l.M12;
                    out.M21 = l.M21;
                }
                m12x *= _b;
                s12x *= _b;
                a12 = sig12 / degree;
                omg12 = lam12 - omg12;
            }
        }

        if ((outmask & Mask.DISTANCE) != 0)
            out.s12 = 0.0 + s12x;
        if ((outmask & Mask.REDUCEDLENGTH) != 0)
            out.m12 = 0.0 + m12x;
        if ((outmask & Mask.AREA) != 0) {
            double // From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
                    salp0 = salp1 * cbet1,
                    calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0
            double alp12;
            if (calp0 != 0 && salp0 != 0) {
                // From Lambda12: tan(bet) = tan(sig) * cos(alp)
                DoublePair sig1 = new DoublePair(sbet1, calp1 * cbet1),
                        sig2 = new DoublePair(sbet2, calp2 * cbet2);
                double k2 = sq(calp0) * _ep2,
                        eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2),
                        // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
                        A4 = sq(_a) * calp0 * salp0 * _e2;
                sinCosNorm(sig1);
                sinCosNorm(sig2);
                double[] C4a = new double[nC4_];
                C4f(eps, C4a);
                double B41 = sinCosSeries(false, sig1.s, sig1.c, C4a, nC4_),
                        B42 = sinCosSeries(false, sig2.s, sig2.c, C4a, nC4_);
                out.S12 = A4 * (B42 - B41);
            } else
                out.S12 = 0;
            if (!meridian &&
                    omg12 < 0.75 * PI &&     // Long difference too big
                    sbet2 - sbet1 < 1.75) {  // Lat difference too big
                // Use tan(Gamma/2) = tan(omg12/2)
                // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
                // with tan(x/2) = sin(x)/(1+cos(x))
                double somg12 = sin(omg12), domg12 = 1 + cos(omg12),
                        dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
                alp12 = 2 * atan2( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
                                   domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) );
            } else {
                // alp12 = alp2 - alp1, used in atan2 so no need to normalize
                double salp12 = salp2 * calp1 - calp2 * salp1,
                        calp12 = calp2 * calp1 + salp2 * salp1;
                // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
                // salp12 = -0 and alp12 = -180.  However this depends on the sign
                // being attached to 0 correctly.  The following ensures the correct
                // behavior.
                if (salp12 == 0 && calp12 < 0) {
                    salp12 = tiny * calp1;
                    calp12 = -1;
                }
                alp12 = atan2(salp12, calp12);
            }
            out.S12 += _c2 * alp12;
            out.S12 *= swapp * lonsign * latsign;
            // convert -0 to +0
            out.S12 += 0;
        }

        // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
        if (swapp < 0) {
            double tmp;
            tmp = salp1; salp1 = salp2; salp2 = tmp;
            tmp = calp1; calp1 = calp2; calp2 = tmp;
            if ((outmask & Mask.GEODESICSCALE) != 0) {
                tmp = out.M12; out.M12 = out.M21; out.M21 = tmp;
            }
        }

        salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
        salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

        if ((outmask & Mask.AZIMUTH) != 0) {
            // minus signs give range [-180, 180). 0- converts -0 to +0.
            out.azi1 = 0 - atan2(-salp1, calp1) / degree;
            out.azi2 = 0 - atan2(-salp2, calp2) / degree;
        }
        // Returned value in [0, 180]
        return a12;
    }

    private static class Lengths {
        double s12b, m12b, m0, M12, M21;
    }

    private Lengths lengths(double eps, double sig12,
                         double ssig1, double csig1, double dn1,
                         double ssig2, double csig2, double dn2,
                         double cbet1, double cbet2,
                         boolean scalep,
                         double[] C1a, double[] C2a) {
        Lengths ret = new Lengths();
        // Return m12b = (reduced length)/_b; also calculate s12b = distance/_b,
        // and m0 = coefficient of secular term in expression for reduced length.
        C1f(eps, C1a);
        C2f(eps, C2a);
        double A1m1 = A1m1f(eps),
                AB1 = (1 + A1m1) * (sinCosSeries(true, ssig2, csig2, C1a, nC1_) -
                                    sinCosSeries(true, ssig1, csig1, C1a, nC1_)),
                A2m1 = A2m1f(eps),
                AB2 = (1 + A2m1) * (sinCosSeries(true, ssig2, csig2, C2a, nC2_) -
                                    sinCosSeries(true, ssig1, csig1, C2a, nC2_));
        ret.m0 = A1m1 - A2m1;
        double J12 = ret.m0 * sig12 + (AB1 - AB2);
        // Missing a factor of _b.
        // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure accurate
        // cancellation in the case of coincident points.
        ret.m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12;
        // Missing a factor of _b
        ret.s12b = (1 + A1m1) * sig12 + AB1;
        if (scalep) {
            double csig12 = csig1 * csig2 + ssig1 * ssig2;
            double t = _ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
            ret.M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
            ret.M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
        }
        return ret;
    }

    public static double astroid(double x, double y) {
        // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
        // This solution is adapted from Geocentric::Reverse.
        double k;
        double p = sq(x),
                q = sq(y),
                r = (p + q - 1) / 6.0;
        if ( !(q == 0 && r <= 0) ) {
            double
                    // Avoid possible division by zero when r = 0 by multiplying equations
                    // for s and t by r^3 and r, resp.
                    S = p * q / 4,            // S = r^3 * s
                    r2 = sq(r),
                    r3 = r * r2,
                    // The discrimant of the quadratic equation for T3.  This is zero on
                    // the evolute curve p^(1/3)+q^(1/3) = 1
                    disc = S * (S + 2 * r3);
            double u = r;
            if (disc >= 0) {
                double T3 = S + r3;
                // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
                // of precision due to cancellation.  The result is unchanged because
                // of the way the T is used in definition of u.
                T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); // T3 = (r * t)^3
                // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
                double T = cbrt(T3); // T = r * t
                // T can be zero; but then r2 / T -> 0.
                u += T + (T != 0 ? r2 / T : 0);
            } else {
                // T is complex, but the way u is defined the result is real.
                double ang = atan2(sqrt(-disc), -(S + r3));
                // There are three possible cube roots.  We choose the root which
                // avoids cancellation.  Note that disc < 0 implies that r < 0.
                u += 2 * r * cos(ang / 3);
            }
            double v = sqrt(sq(u) + q),    // guaranteed positive
                    // Avoid loss of accuracy when u < 0.
                    uv = u < 0 ? q / (v - u) : u + v, // u+v, guaranteed positive
                    w = (uv - q) / (2 * v);           // positive?
            // Rearrange expression for k to avoid loss of accuracy due to
            // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
            k = uv / (sqrt(uv + sq(w)) + w);   // guaranteed positive
        } else {               // q == 0 && r <= 0
            // y = 0 with |x| <= 1.  Handle this case directly.
            // for y small, positive root is k = abs(y)/sqrt(1-x^2)
            k = 0;
        }
        return k;
    }

    private double inverseStart(double sbet1, double cbet1, double dn1,
                                    double sbet2, double cbet2, double dn2,
                                    double lam12,
                                    DoublePair alp1,
                                    // Only updated if return val >= 0
                                    DoublePair alp2,
                                    // Scratch areas of the right size
                                    double C1a[], double C2a[]) {
        // Return a starting point for Newton's method in salp1 and calp1 (function
        // value is -1).  If Newton's method doesn't need to be used, return also
        // salp2 and calp2 and function value is sig12.
        double sig12 = -1,               // Return value
                // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
                sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
                cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
        double sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
        boolean shortline = cbet12 >= 0 && sbet12 < 0.5 &&
                lam12 <= PI / 6;
        double omg12 = !shortline ? lam12 : lam12 / (_f1 * (dn1 + dn2) / 2),
                somg12 = sin(omg12), comg12 = cos(omg12);

        alp1.s = cbet2 * somg12;
        alp1.c = comg12 >= 0 ?
                sbet12 + cbet2 * sbet1 * sq(somg12) / (1 + comg12) :
                sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
        double ssig12 = hypot(alp1.s, alp1.c),
                csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;
        if (shortline && ssig12 < _etol2) {
            // really short lines
            alp2.s = cbet1 * somg12;
            alp2.c = sbet12 - cbet1 * sbet2 * sq(somg12) / (1 + comg12);
            sinCosNorm(alp2);
            // Set return value
            sig12 = atan2(ssig12, csig12);
        } else if (abs(_n) > 0.1 || // Skip astroid calc if too eccentric
                csig12 >= 0 ||
                ssig12 >= 6 * abs(_n) * PI * sq(cbet1)) {
            // Nothing to do, zeroth order spherical approximation is OK
        } else {
            // Scale lam12 and bet2 to x, y coordinate system where antipodal point
            // is at origin and singular point is at y = 0, x = -1.
            double y, lamscale, betscale;
            // Volatile declaration needed to fix inverse case
            // 56.320923501171 0 -56.320923501171 179.664747671772880215
            // which otherwise fails with g++ 4.4.4 x86 -O3
            double x;
            if (_f >= 0) {            // In fact f == 0 does not get here
                // x = dlong, y = dlat
                {
                    double k2 = sq(sbet1) * _ep2,
                            eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
                    lamscale = _f * cbet1 * A3f(eps) * PI;
                }
                betscale = lamscale * cbet1;
                x = (lam12 - PI) / lamscale;
                y = sbet12a / betscale;
            } else {                  // _f < 0
                // x = dlat, y = dlong
                double cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
                        bet12a = atan2(sbet12a, cbet12a);
                double m12b, m0;
                // In the case of lon12 = 180, this repeats a calculation made in
                // Inverse.
                Lengths l = lengths(_n, PI + bet12a,
                        sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                        cbet1, cbet2, false,
                        C1a, C2a);
                x = -1 + l.m12b / (cbet1 * cbet2 * l.m0 * PI);
                betscale = x < -0.01 ? sbet12a / x :
                        -_f * sq(cbet1) * PI;
                lamscale = betscale / cbet1;
                y = (lam12 - PI) / lamscale;
            }

            if (y > -tol1 && x > -1 - xthresh) {
                // strip near cut
                // Need real(x) here to cast away the volatility of x for min/max
                if (_f >= 0) {
                    alp1.s = min(1.0, -x); alp1.c = -sqrt(1 - sq(alp1.s));
                } else {
                    alp1.c = max(x > -tol1 ? 0.0 : -1.0, x);
                    alp1.s = sqrt(1 - alp1.c);
                }
            } else {
                // Estimate alp1, by solving the astroid problem.
                //
                // Could estimate alpha1 = theta + pi/2, directly, i.e.,
                //   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
                //   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
                //
                // However, it's better to estimate omg12 from astroid and use
                // spherical formula to compute alp1.  This reduces the mean number of
                // Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
                // (min 0 max 5).  The changes in the number of iterations are as
                // follows:
                //
                // change percent
                //    1       5
                //    0      78
                //   -1      16
                //   -2       0.6
                //   -3       0.04
                //   -4       0.002
                //
                // The histogram of iterations is (m = number of iterations estimating
                // alp1 directly, n = number of iterations estimating via omg12, total
                // number of trials = 148605):
                //
                //  iter    m      n
                //    0   148    186
                //    1 13046  13845
                //    2 93315 102225
                //    3 36189  32341
                //    4  5396      7
                //    5   455      1
                //    6    56      0
                //
                // Because omg12 is near pi, estimate work with omg12a = pi - omg12
                double k = astroid(x, y);
                double omg12a = lamscale * ( _f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k );
                somg12 = sin(omg12a); comg12 = -cos(omg12a);
                // Update spherical estimate of alp1 using omg12 instead of lam12
                alp1.s = cbet2 * somg12;
                alp1.c = sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
            }
        }
        if (alp1.s > 0) // Sanity check on starting guess
            sinCosNorm(alp1);
        else {
            alp1.s = 1;
            alp1.c = 0;
        }
        return sig12;
    }

    private static class LambdaResult {
        public double salp2, calp2, sig12, eps, domg12, dlam12;
        DoublePair sig1, sig2;
        public LambdaResult() { sig1 = new DoublePair(0.0, 0.0); sig2 = new DoublePair(0.0, 0.0); }
    }

    private double lambda12(double sbet1, double cbet1, double dn1,
                            double sbet2, double cbet2, double dn2,
                            double salp1, double calp1,
                            boolean diffp, LambdaResult out,
                            // Scratch areas of the right size
                            double[] C1a, double[] C2a, double[] C3a) {
        if (sbet1 == 0 && calp1 == 0) {
            // Break degeneracy of equatorial line.  This case has already been
            // handled.
            calp1 = -tiny;
        }

        double // sin(alp1) * cos(bet1) = sin(alp0)
                salp0 = salp1 * cbet1,
                calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0

        double somg1, comg1, somg2, comg2, omg12, lam12;
        // tan(bet1) = tan(sig1) * cos(alp1)
        // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
        out.sig1.s = sbet1; somg1 = salp0 * sbet1;
        out.sig1.c = comg1 = calp1 * cbet1;
        sinCosNorm(out.sig1);
        // SinCosNorm(somg1, comg1); -- don't need to normalize!

        // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
        // about this case, since this can yield singularities in the Newton
        // iteration.
        // sin(alp2) * cos(bet2) = sin(alp0)
        out.salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
        // calp2 = sqrt(1 - sq(salp2))
        //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
        // and subst for calp0 and rearrange to give (choose positive sqrt
        // to give alp2 in [0, pi/2]).
        out.calp2 = cbet2 != cbet1 || abs(sbet2) != -sbet1 ?
                sqrt(sq(calp1 * cbet1) +
                        (cbet1 < -sbet1 ?
                                (cbet2 - cbet1) * (cbet1 + cbet2) :
                                (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
                abs(calp1);
        // tan(bet2) = tan(sig2) * cos(alp2)
        // tan(omg2) = sin(alp0) * tan(sig2).
        out.sig2.s = sbet2; somg2 = salp0 * sbet2;
        out.sig2.c = comg2 = out.calp2 * cbet2;
        sinCosNorm(out.sig2);
        // SinCosNorm(somg2, comg2); -- don't need to normalize!

        // sig12 = sig2 - sig1, limit to [0, pi]
        out.sig12 = atan2(max(out.sig1.c * out.sig2.s - out.sig1.s * out.sig2.c, 0.0),
                out.sig1.c * out.sig2.c + out.sig1.s * out.sig2.s);

        // omg12 = omg2 - omg1, limit to [0, pi]
        omg12 = atan2(max(comg1 * somg2 - somg1 * comg2, 0.0),
                comg1 * comg2 + somg1 * somg2);
        double B312, h0;
        double k2 = sq(calp0) * _ep2;
        out.eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
        C3f(out.eps, C3a);
        B312 = (sinCosSeries(true, out.sig2.s, out.sig2.c, C3a, nC3_-1) -
                sinCosSeries(true, out.sig1.s, out.sig1.c, C3a, nC3_-1));
        h0 = -_f * A3f(out.eps);
        out.domg12 = salp0 * h0 * (out.sig12 + B312);
        lam12 = omg12 + out.domg12;
        if (diffp) {
            if (out.calp2 == 0)
                out.dlam12 = - 2 * _f1 * dn1 / sbet1;
            else {
                Lengths l = lengths(out.eps, out.sig12, out.sig1.s, out.sig1.c, dn1, out.sig2.s, out.sig2.c, dn2,
                        cbet1, cbet2, false, C1a, C2a);
                out.dlam12 = l.m12b * _f1 / (out.calp2 * cbet2);
            }
        }

        return lam12;
    }

    double A3f(double eps) {
        // Evaluate sum(_A3x[k] * eps^k, k, 0, nA3x_-1) by Horner's method
        double v = 0;
        for (int i = nA3x_; i != 0; )
            v = eps * v + _A3x[--i];
        return v;
    }

    void C3f(double eps, double[] c) {
        // Evaluate C3 coeffs by Horner's method
        // Elements c[1] thru c[nC3_ - 1] are set
        for (int j = nC3x_, k = nC3_ - 1; k != 0; ) {
            double t = 0;
            for (int i = nC3_ - k; i != 0; --i)
                t = eps * t + _C3x[--j];
            c[k--] = t;
        }

        double mult = 1;
        for (int k = 1; k < nC3_; ) {
            mult *= eps;
            c[k++] *= mult;
        }
    }

    void C4f(double eps, double c[]) {
        // Evaluate C4 coeffs by Horner's method
        // Elements c[0] thru c[nC4_ - 1] are set
        for (int j = nC4x_, k = nC4_; k != 0; ) {
            double t = 0;
            for (int i = nC4_ - k + 1; i != 0; --i)
                t = eps * t + _C4x[--j];
            c[--k] = t;
        }

        double mult = 1;
        for (int k = 1; k < nC4_; ) {
            mult *= eps;
            c[k++] *= mult;
        }
    }

    static double A1m1f(double eps) {
        final double eps2 = sq(eps), t;
        switch (nA1_/2) {
        case 0:
            t = 0;
            break;
        case 1:
            t = eps2/4;
            break;
        case 2:
            t = eps2*(eps2+16)/64;
            break;
        case 3:
            t = eps2*(eps2*(eps2+4)+64)/256;
            break;
        case 4:
            t = eps2*(eps2*(eps2*(25*eps2+64)+256)+4096)/16384;
            break;
        default:
            throw new IllegalStateException("Bad value of nA1_");
        }
        return (t + eps) / (1 - eps);
    }

    static void C1f(double eps, double[] c) {
        final double eps2 = sq(eps);
        double d = eps;
        switch (nC1_) {
        case 0:
            break;
        case 1:
            c[1] = -d/2;
            break;
        case 2:
            c[1] = -d/2;
            d *= eps;
            c[2] = -d/16;
            break;
        case 3:
            c[1] = d*(3*eps2-8)/16;
            d *= eps;
            c[2] = -d/16;
            d *= eps;
            c[3] = -d/48;
            break;
        case 4:
            c[1] = d*(3*eps2-8)/16;
            d *= eps;
            c[2] = d*(eps2-2)/32;
            d *= eps;
            c[3] = -d/48;
            d *= eps;
            c[4] = -5*d/512;
            break;
        case 5:
            c[1] = d*((6-eps2)*eps2-16)/32;
            d *= eps;
            c[2] = d*(eps2-2)/32;
            d *= eps;
            c[3] = d*(9*eps2-16)/768;
            d *= eps;
            c[4] = -5*d/512;
            d *= eps;
            c[5] = -7*d/1280;
            break;
        case 6:
            c[1] = d*((6-eps2)*eps2-16)/32;
            d *= eps;
            c[2] = d*((64-9*eps2)*eps2-128)/2048;
            d *= eps;
            c[3] = d*(9*eps2-16)/768;
            d *= eps;
            c[4] = d*(3*eps2-5)/512;
            d *= eps;
            c[5] = -7*d/1280;
            d *= eps;
            c[6] = -7*d/2048;
            break;
        case 7:
            c[1] = d*(eps2*(eps2*(19*eps2-64)+384)-1024)/2048;
            d *= eps;
            c[2] = d*((64-9*eps2)*eps2-128)/2048;
            d *= eps;
            c[3] = d*((72-9*eps2)*eps2-128)/6144;
            d *= eps;
            c[4] = d*(3*eps2-5)/512;
            d *= eps;
            c[5] = d*(35*eps2-56)/10240;
            d *= eps;
            c[6] = -7*d/2048;
            d *= eps;
            c[7] = -33*d/14336;
            break;
        case 8:
            c[1] = d*(eps2*(eps2*(19*eps2-64)+384)-1024)/2048;
            d *= eps;
            c[2] = d*(eps2*(eps2*(7*eps2-18)+128)-256)/4096;
            d *= eps;
            c[3] = d*((72-9*eps2)*eps2-128)/6144;
            d *= eps;
            c[4] = d*((96-11*eps2)*eps2-160)/16384;
            d *= eps;
            c[5] = d*(35*eps2-56)/10240;
            d *= eps;
            c[6] = d*(9*eps2-14)/4096;
            d *= eps;
            c[7] = -33*d/14336;
            d *= eps;
            c[8] = -429*d/262144;
            break;
        default:
            throw new IllegalStateException("Bad value of nC1_");
        }
    }

    // The coefficients C1p[l] in the Fourier expansion of B1p
    static void C1pf(double eps, double[] c) {
        double eps2 = sq(eps),
                d = eps;
        switch (nC1p_) {
        case 0:
            break;
        case 1:
            c[1] = d/2;
            break;
        case 2:
            c[1] = d/2;
            d *= eps;
            c[2] = 5*d/16;
            break;
        case 3:
            c[1] = d*(16-9*eps2)/32;
            d *= eps;
            c[2] = 5*d/16;
            d *= eps;
            c[3] = 29*d/96;
            break;
        case 4:
            c[1] = d*(16-9*eps2)/32;
            d *= eps;
            c[2] = d*(30-37*eps2)/96;
            d *= eps;
            c[3] = 29*d/96;
            d *= eps;
            c[4] = 539*d/1536;
            break;
        case 5:
            c[1] = d*(eps2*(205*eps2-432)+768)/1536;
            d *= eps;
            c[2] = d*(30-37*eps2)/96;
            d *= eps;
            c[3] = d*(116-225*eps2)/384;
            d *= eps;
            c[4] = 539*d/1536;
            d *= eps;
            c[5] = 3467*d/7680;
            break;
        case 6:
            c[1] = d*(eps2*(205*eps2-432)+768)/1536;
            d *= eps;
            c[2] = d*(eps2*(4005*eps2-4736)+3840)/12288;
            d *= eps;
            c[3] = d*(116-225*eps2)/384;
            d *= eps;
            c[4] = d*(2695-7173*eps2)/7680;
            d *= eps;
            c[5] = 3467*d/7680;
            d *= eps;
            c[6] = 38081*d/61440;
            break;
        case 7:
            c[1] = d*(eps2*((9840-4879*eps2)*eps2-20736)+36864)/73728;
            d *= eps;
            c[2] = d*(eps2*(4005*eps2-4736)+3840)/12288;
            d *= eps;
            c[3] = d*(eps2*(8703*eps2-7200)+3712)/12288;
            d *= eps;
            c[4] = d*(2695-7173*eps2)/7680;
            d *= eps;
            c[5] = d*(41604-141115*eps2)/92160;
            d *= eps;
            c[6] = 38081*d/61440;
            d *= eps;
            c[7] = 459485*d/516096;
            break;
        case 8:
            c[1] = d*(eps2*((9840-4879*eps2)*eps2-20736)+36864)/73728;
            d *= eps;
            c[2] = d*(eps2*((120150-86171*eps2)*eps2-142080)+115200)/368640;
            d *= eps;
            c[3] = d*(eps2*(8703*eps2-7200)+3712)/12288;
            d *= eps;
            c[4] = d*(eps2*(1082857*eps2-688608)+258720)/737280;
            d *= eps;
            c[5] = d*(41604-141115*eps2)/92160;
            d *= eps;
            c[6] = d*(533134-2200311*eps2)/860160;
            d *= eps;
            c[7] = 459485*d/516096;
            d *= eps;
            c[8] = 109167851*d/82575360;
            break;
        default:
            throw new IllegalStateException("Bad value of nC1p_");
        }
    }

    // The scale factor A2-1 = mean value of (d/dsigma)I2 - 1
    static double A2m1f(double eps) {
        double eps2 = sq(eps),
                t;
        switch (nA2_/2) {
        case 0:
            t = 0;
            break;
        case 1:
            t = eps2/4;
            break;
        case 2:
            t = eps2*(9*eps2+16)/64;
            break;
        case 3:
            t = eps2*(eps2*(25*eps2+36)+64)/256;
            break;
        case 4:
            t = eps2*(eps2*(eps2*(1225*eps2+1600)+2304)+4096)/16384;
            break;
        default:
            throw new IllegalStateException("Bad value of nA2_");
        }
        return t * (1 - eps) - eps;
    }

    // The coefficients C2[l] in the Fourier expansion of B2
    static void C2f(double eps, double[] c) {
        double
                eps2 = sq(eps),
                d = eps;
        switch (nC2_) {
        case 0:
            break;
        case 1:
            c[1] = d/2;
            break;
        case 2:
            c[1] = d/2;
            d *= eps;
            c[2] = 3*d/16;
            break;
        case 3:
            c[1] = d*(eps2+8)/16;
            d *= eps;
            c[2] = 3*d/16;
            d *= eps;
            c[3] = 5*d/48;
            break;
        case 4:
            c[1] = d*(eps2+8)/16;
            d *= eps;
            c[2] = d*(eps2+6)/32;
            d *= eps;
            c[3] = 5*d/48;
            d *= eps;
            c[4] = 35*d/512;
            break;
        case 5:
            c[1] = d*(eps2*(eps2+2)+16)/32;
            d *= eps;
            c[2] = d*(eps2+6)/32;
            d *= eps;
            c[3] = d*(15*eps2+80)/768;
            d *= eps;
            c[4] = 35*d/512;
            d *= eps;
            c[5] = 63*d/1280;
            break;
        case 6:
            c[1] = d*(eps2*(eps2+2)+16)/32;
            d *= eps;
            c[2] = d*(eps2*(35*eps2+64)+384)/2048;
            d *= eps;
            c[3] = d*(15*eps2+80)/768;
            d *= eps;
            c[4] = d*(7*eps2+35)/512;
            d *= eps;
            c[5] = 63*d/1280;
            d *= eps;
            c[6] = 77*d/2048;
            break;
        case 7:
            c[1] = d*(eps2*(eps2*(41*eps2+64)+128)+1024)/2048;
            d *= eps;
            c[2] = d*(eps2*(35*eps2+64)+384)/2048;
            d *= eps;
            c[3] = d*(eps2*(69*eps2+120)+640)/6144;
            d *= eps;
            c[4] = d*(7*eps2+35)/512;
            d *= eps;
            c[5] = d*(105*eps2+504)/10240;
            d *= eps;
            c[6] = 77*d/2048;
            d *= eps;
            c[7] = 429*d/14336;
            break;
        case 8:
            c[1] = d*(eps2*(eps2*(41*eps2+64)+128)+1024)/2048;
            d *= eps;
            c[2] = d*(eps2*(eps2*(47*eps2+70)+128)+768)/4096;
            d *= eps;
            c[3] = d*(eps2*(69*eps2+120)+640)/6144;
            d *= eps;
            c[4] = d*(eps2*(133*eps2+224)+1120)/16384;
            d *= eps;
            c[5] = d*(105*eps2+504)/10240;
            d *= eps;
            c[6] = d*(33*eps2+154)/4096;
            d *= eps;
            c[7] = 429*d/14336;
            d *= eps;
            c[8] = 6435*d/262144;
            break;
        default:
            throw new IllegalStateException("Bad value of nC2_");
        }
    }

    private void A3coeff() {
        switch (nA3_) {
        case 0:
            break;
        case 1:
            _A3x[0] = 1;
            break;
        case 2:
            _A3x[0] = 1;
            _A3x[1] = -1/2.0;
            break;
        case 3:
            _A3x[0] = 1;
            _A3x[1] = (_n-1)/2;
            _A3x[2] = -1/4.0;
            break;
        case 4:
            _A3x[0] = 1;
            _A3x[1] = (_n-1)/2;
            _A3x[2] = (-_n-2)/8;
            _A3x[3] = -1/16.0;
            break;
        case 5:
            _A3x[1] = (_n-1)/2;
            _A3x[2] = (_n*(3*_n-1)-2)/8;
            _A3x[3] = (-3*_n-1)/16;
            _A3x[4] = -3/64.0;
            break;
        case 6:
            _A3x[0] = 1;
            _A3x[1] = (_n-1)/2;
            _A3x[2] = (_n*(3*_n-1)-2)/8;
            _A3x[3] = ((-_n-3)*_n-1)/16;
            _A3x[4] = (-2*_n-3)/64;
            _A3x[5] = -3/128.0;
            break;
        case 7:
            _A3x[0] = 1;
            _A3x[1] = (_n-1)/2;
            _A3x[2] = (_n*(3*_n-1)-2)/8;
            _A3x[3] = (_n*(_n*(5*_n-1)-3)-1)/16;
            _A3x[4] = ((-10*_n-2)*_n-3)/64;
            _A3x[5] = (-5*_n-3)/128;
            _A3x[6] = -5/256.0;
            break;
        case 8:
            _A3x[0] = 1;
            _A3x[1] = (_n-1)/2;
            _A3x[2] = (_n*(3*_n-1)-2)/8;
            _A3x[3] = (_n*(_n*(5*_n-1)-3)-1)/16;
            _A3x[4] = (_n*((-5*_n-20)*_n-4)-6)/128;
            _A3x[5] = ((-5*_n-10)*_n-6)/256;
            _A3x[6] = (-15*_n-20)/1024;
            _A3x[7] = -25/2048.0;
            break;
        default:
            throw new IllegalStateException("Bad value of nA3_");

        }
    }

    private void C3coeff() {
        switch (nC3_) {
        case 0:
            break;
        case 1:
            break;
        case 2:
            _C3x[0] = 1/4.0;
            break;
        case 3:
            _C3x[0] = (1-_n)/4;
            _C3x[1] = 1/8.0;
            _C3x[2] = 1/16.0;
            break;
        case 4:
            _C3x[0] = (1-_n)/4;
            _C3x[1] = 1/8.0;
            _C3x[2] = 3/64.0;
            _C3x[3] = (2-3*_n)/32;
            _C3x[4] = 3/64.0;
            _C3x[5] = 5/192.0;
            break;
        case 5:
            _C3x[0] = (1-_n)/4;
            _C3x[1] = (1-_n*_n)/8;
            _C3x[2] = (3*_n+3)/64;
            _C3x[3] = 5/128.0;
            _C3x[4] = ((_n-3)*_n+2)/32;
            _C3x[5] = (3-2*_n)/64;
            _C3x[6] = 3/128.0;
            _C3x[7] = (5-9*_n)/192;
            _C3x[8] = 3/128.0;
            _C3x[9] = 7/512.0;
            break;
        case 6:
            _C3x[0] = (1-_n)/4;
            _C3x[1] = (1-_n*_n)/8;
            _C3x[2] = ((3-_n)*_n+3)/64;
            _C3x[3] = (2*_n+5)/128;
            _C3x[4] = 3/128.0;
            _C3x[5] = ((_n-3)*_n+2)/32;
            _C3x[6] = ((-3*_n-2)*_n+3)/64;
            _C3x[7] = (_n+3)/128;
            _C3x[8] = 5/256.0;
            _C3x[9] = (_n*(5*_n-9)+5)/192;
            _C3x[10] = (9-10*_n)/384;
            _C3x[11] = 7/512.0;
            _C3x[12] = (7-14*_n)/512;
            _C3x[13] = 7/512.0;
            _C3x[14] = 21/2560.0;
            break;
        case 7:
            _C3x[0] = (1-_n)/4;
            _C3x[1] = (1-_n*_n)/8;
            _C3x[2] = (_n*((-5*_n-1)*_n+3)+3)/64;
            _C3x[3] = (_n*(2*_n+2)+5)/128;
            _C3x[4] = (11*_n+12)/512;
            _C3x[5] = 21/1024.0;
            _C3x[6] = ((_n-3)*_n+2)/32;
            _C3x[7] = (_n*(_n*(2*_n-3)-2)+3)/64;
            _C3x[8] = ((2-9*_n)*_n+6)/256;
            _C3x[9] = (_n+5)/256;
            _C3x[10] = 27/2048.0;
            _C3x[11] = (_n*((5-_n)*_n-9)+5)/192;
            _C3x[12] = ((-6*_n-10)*_n+9)/384;
            _C3x[13] = (21-4*_n)/1536;
            _C3x[14] = 3/256.0;
            _C3x[15] = (_n*(10*_n-14)+7)/512;
            _C3x[16] = (7-10*_n)/512;
            _C3x[17] = 9/1024.0;
            _C3x[18] = (21-45*_n)/2560;
            _C3x[19] = 9/1024.0;
            _C3x[20] = 11/2048.0;
            break;
        case 8:
            _C3x[0] = (1-_n)/4;
            _C3x[1] = (1-_n*_n)/8;
            _C3x[2] = (_n*((-5*_n-1)*_n+3)+3)/64;
            _C3x[3] = (_n*((2-2*_n)*_n+2)+5)/128;
            _C3x[4] = (_n*(3*_n+11)+12)/512;
            _C3x[5] = (10*_n+21)/1024;
            _C3x[6] = 243/16384.0;
            _C3x[7] = ((_n-3)*_n+2)/32;
            _C3x[8] = (_n*(_n*(2*_n-3)-2)+3)/64;
            _C3x[9] = (_n*((-6*_n-9)*_n+2)+6)/256;
            _C3x[10] = ((1-2*_n)*_n+5)/256;
            _C3x[11] = (69*_n+108)/8192;
            _C3x[12] = 187/16384.0;
            _C3x[13] = (_n*((5-_n)*_n-9)+5)/192;
            _C3x[14] = (_n*(_n*(10*_n-6)-10)+9)/384;
            _C3x[15] = ((-77*_n-8)*_n+42)/3072;
            _C3x[16] = (12-_n)/1024;
            _C3x[17] = 139/16384.0;
            _C3x[18] = (_n*((20-7*_n)*_n-28)+14)/1024;
            _C3x[19] = ((-7*_n-40)*_n+28)/2048;
            _C3x[20] = (72-43*_n)/8192;
            _C3x[21] = 127/16384.0;
            _C3x[22] = (_n*(75*_n-90)+42)/5120;
            _C3x[23] = (9-15*_n)/1024;
            _C3x[24] = 99/16384.0;
            _C3x[25] = (44-99*_n)/8192;
            _C3x[26] = 99/16384.0;
            _C3x[27] = 429/114688.0;
            break;
        default:
            throw new IllegalStateException("Bad value of nC3_");
        }
    }

    private void C4coeff() {
        switch (nC4_) {
        case 0:
            break;
        case 1:
            _C4x[0] = 2/3.0;
            break;
        case 2:
            _C4x[0] = (10-4*_n)/15;
            _C4x[1] = -1/5.0;
            _C4x[2] = 1/45.0;
            break;
        case 3:
            _C4x[0] = (_n*(8*_n-28)+70)/105;
            _C4x[1] = (16*_n-7)/35;
            _C4x[2] = -2/105.0;
            _C4x[3] = (7-16*_n)/315;
            _C4x[4] = -2/105.0;
            _C4x[5] = 4/525.0;
            break;
        case 4:
            _C4x[0] = (_n*(_n*(4*_n+24)-84)+210)/315;
            _C4x[1] = ((48-32*_n)*_n-21)/105;
            _C4x[2] = (-32*_n-6)/315;
            _C4x[3] = 11/315.0;
            _C4x[4] = (_n*(32*_n-48)+21)/945;
            _C4x[5] = (64*_n-18)/945;
            _C4x[6] = -1/105.0;
            _C4x[7] = (12-32*_n)/1575;
            _C4x[8] = -8/1575.0;
            _C4x[9] = 8/2205.0;
            break;
        case 5:
            _C4x[0] = (_n*(_n*(_n*(16*_n+44)+264)-924)+2310)/3465;
            _C4x[1] = (_n*(_n*(48*_n-352)+528)-231)/1155;
            _C4x[2] = (_n*(1088*_n-352)-66)/3465;
            _C4x[3] = (121-368*_n)/3465;
            _C4x[4] = 4/1155.0;
            _C4x[5] = (_n*((352-48*_n)*_n-528)+231)/10395;
            _C4x[6] = ((704-896*_n)*_n-198)/10395;
            _C4x[7] = (80*_n-99)/10395;
            _C4x[8] = 4/1155.0;
            _C4x[9] = (_n*(320*_n-352)+132)/17325;
            _C4x[10] = (384*_n-88)/17325;
            _C4x[11] = -8/1925.0;
            _C4x[12] = (88-256*_n)/24255;
            _C4x[13] = -16/8085.0;
            _C4x[14] = 64/31185.0;
            break;
        case 6:
            _C4x[0] = (_n*(_n*(_n*(_n*(100*_n+208)+572)+3432)-12012)+30030)/45045;
            _C4x[1] = (_n*(_n*(_n*(64*_n+624)-4576)+6864)-3003)/15015;
            _C4x[2] = (_n*((14144-10656*_n)*_n-4576)-858)/45045;
            _C4x[3] = ((-224*_n-4784)*_n+1573)/45045;
            _C4x[4] = (1088*_n+156)/45045;
            _C4x[5] = 97/15015.0;
            _C4x[6] = (_n*(_n*((-64*_n-624)*_n+4576)-6864)+3003)/135135;
            _C4x[7] = (_n*(_n*(5952*_n-11648)+9152)-2574)/135135;
            _C4x[8] = (_n*(5792*_n+1040)-1287)/135135;
            _C4x[9] = (468-2944*_n)/135135;
            _C4x[10] = 1/9009.0;
            _C4x[11] = (_n*((4160-1440*_n)*_n-4576)+1716)/225225;
            _C4x[12] = ((4992-8448*_n)*_n-1144)/225225;
            _C4x[13] = (1856*_n-936)/225225;
            _C4x[14] = 8/10725.0;
            _C4x[15] = (_n*(3584*_n-3328)+1144)/315315;
            _C4x[16] = (1024*_n-208)/105105;
            _C4x[17] = -136/63063.0;
            _C4x[18] = (832-2560*_n)/405405;
            _C4x[19] = -128/135135.0;
            _C4x[20] = 128/99099.0;
            break;
        case 7:
            _C4x[0] = (_n*(_n*(_n*(_n*(_n*(56*_n+100)+208)+572)+3432)-12012)+30030)/
                    45045;
            _C4x[1] = (_n*(_n*(_n*(_n*(16*_n+64)+624)-4576)+6864)-3003)/15015;
            _C4x[2] = (_n*(_n*(_n*(1664*_n-10656)+14144)-4576)-858)/45045;
            _C4x[3] = (_n*(_n*(10736*_n-224)-4784)+1573)/45045;
            _C4x[4] = ((1088-4480*_n)*_n+156)/45045;
            _C4x[5] = (291-464*_n)/45045;
            _C4x[6] = 10/9009.0;
            _C4x[7] = (_n*(_n*(_n*((-16*_n-64)*_n-624)+4576)-6864)+3003)/135135;
            _C4x[8] = (_n*(_n*((5952-768*_n)*_n-11648)+9152)-2574)/135135;
            _C4x[9] = (_n*((5792-10704*_n)*_n+1040)-1287)/135135;
            _C4x[10] = (_n*(3840*_n-2944)+468)/135135;
            _C4x[11] = (112*_n+15)/135135;
            _C4x[12] = 10/9009.0;
            _C4x[13] = (_n*(_n*(_n*(128*_n-1440)+4160)-4576)+1716)/225225;
            _C4x[14] = (_n*(_n*(6784*_n-8448)+4992)-1144)/225225;
            _C4x[15] = (_n*(1664*_n+1856)-936)/225225;
            _C4x[16] = (168-1664*_n)/225225;
            _C4x[17] = -4/25025.0;
            _C4x[18] = (_n*((3584-1792*_n)*_n-3328)+1144)/315315;
            _C4x[19] = ((1024-2048*_n)*_n-208)/105105;
            _C4x[20] = (1792*_n-680)/315315;
            _C4x[21] = 64/315315.0;
            _C4x[22] = (_n*(3072*_n-2560)+832)/405405;
            _C4x[23] = (2048*_n-384)/405405;
            _C4x[24] = -512/405405.0;
            _C4x[25] = (640-2048*_n)/495495;
            _C4x[26] = -256/495495.0;
            _C4x[27] = 512/585585.0;
            break;
        case 8:
            _C4x[0] = (_n*(_n*(_n*(_n*(_n*(_n*(588*_n+952)+1700)+3536)+9724)+58344)-
                    204204)+510510)/765765;
            _C4x[1] = (_n*(_n*(_n*(_n*(_n*(96*_n+272)+1088)+10608)-77792)+116688)-
                    51051)/255255;
            _C4x[2] = (_n*(_n*(_n*(_n*(3232*_n+28288)-181152)+240448)-77792)-14586)/
                    765765;
            _C4x[3] = (_n*(_n*((182512-154048*_n)*_n-3808)-81328)+26741)/765765;
            _C4x[4] = (_n*(_n*(12480*_n-76160)+18496)+2652)/765765;
            _C4x[5] = (_n*(20960*_n-7888)+4947)/765765;
            _C4x[6] = (4192*_n+850)/765765;
            _C4x[7] = 193/85085.0;
            _C4x[8] = (_n*(_n*(_n*(_n*((-96*_n-272)*_n-1088)-10608)+77792)-116688)+
                    51051)/2297295;
            _C4x[9] = (_n*(_n*(_n*((-1344*_n-13056)*_n+101184)-198016)+155584)-43758)/
                    2297295;
            _C4x[10] = (_n*(_n*(_n*(103744*_n-181968)+98464)+17680)-21879)/2297295;
            _C4x[11] = (_n*(_n*(52608*_n+65280)-50048)+7956)/2297295;
            _C4x[12] = ((1904-39840*_n)*_n+255)/2297295;
            _C4x[13] = (510-1472*_n)/459459;
            _C4x[14] = 349/2297295.0;
            _C4x[15] = (_n*(_n*(_n*(_n*(160*_n+2176)-24480)+70720)-77792)+29172)/
                    3828825;
            _C4x[16] = (_n*(_n*((115328-41472*_n)*_n-143616)+84864)-19448)/3828825;
            _C4x[17] = (_n*((28288-126528*_n)*_n+31552)-15912)/3828825;
            _C4x[18] = (_n*(64256*_n-28288)+2856)/3828825;
            _C4x[19] = (-928*_n-612)/3828825;
            _C4x[20] = 464/1276275.0;
            _C4x[21] = (_n*(_n*(_n*(7168*_n-30464)+60928)-56576)+19448)/5360355;
            _C4x[22] = (_n*(_n*(35840*_n-34816)+17408)-3536)/1786785;
            _C4x[23] = ((30464-2560*_n)*_n-11560)/5360355;
            _C4x[24] = (1088-16384*_n)/5360355;
            _C4x[25] = -16/97461.0;
            _C4x[26] = (_n*((52224-32256*_n)*_n-43520)+14144)/6891885;
            _C4x[27] = ((34816-77824*_n)*_n-6528)/6891885;
            _C4x[28] = (26624*_n-8704)/6891885;
            _C4x[29] = 128/2297295.0;
            _C4x[30] = (_n*(45056*_n-34816)+10880)/8423415;
            _C4x[31] = (24576*_n-4352)/8423415;
            _C4x[32] = -6784/8423415.0;
            _C4x[33] = (8704-28672*_n)/9954945;
            _C4x[34] = -1024/3318315.0;
            _C4x[35] = 1024/1640925.0;
            break;
        default:
            throw new IllegalStateException("Bad value of nC4_");
        }
    }
}