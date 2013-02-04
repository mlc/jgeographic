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
import static com.github.mlc.jgeographic.Geodesic.Mask.*;

public final class GeodesicLine {
    static final int nC1_ = Geodesic.nC1_;
    static final int nC1p_ = Geodesic.nC1p_;
    static final int nC2_ = Geodesic.nC2_;
    static final int nC3_ = Geodesic.nC3_;
    static final int nC4_ = Geodesic.nC4_;

    final double _lat1, _lon1, _azi1;
    final double  _a, _f, _b, _c2, _f1, _salp0, _calp0, _k2,
            _salp1, _calp1, _ssig1, _csig1, _dn1;
    double _stau1, _ctau1, _somg1, _comg1,
            _A1m1, _A2m1, _A3c, _B11, _B21, _B31, _A4, _B41;
    // index zero elements of _C1a, _C1pa, _C2a, _C3a are unused
    final double[] _C1a = new double[nC1_ + 1], _C1pa = new double[nC1p_ + 1],
            _C2a = new double[nC2_ + 1], _C3a = new double[nC3_],
            _C4a = new double[nC4_];
    final int _caps;

    static final int CAP_NONE = Geodesic.CAP_NONE,
            CAP_C1   = Geodesic.CAP_C1,
            CAP_C1p  = Geodesic.CAP_C1p,
            CAP_C2   = Geodesic.CAP_C2,
            CAP_C3   = Geodesic.CAP_C3,
            CAP_C4   = Geodesic.CAP_C4,
            CAP_ALL  = Geodesic.CAP_ALL,
            OUT_ALL  = Geodesic.OUT_ALL;

    public GeodesicLine(Geodesic g, double lat1, double lon1, double azi1, int caps) {
        _a = g._a;
        _f = g._f;
        _b = g._b;
        _c2 = g._c2;
        _f1 = g._f1;
        _caps = caps | LATITUDE | LONGITUDE;
        azi1 = Geodesic.angRound(azi1);
        lon1 = angNormalize(lon1);
        _lat1 = lat1;
        _lon1 = lon1;
        _azi1 = azi1;
        double alp1 = azi1 * degree;
        _salp1 =     azi1  == -180 ? 0 : sin(alp1);
        _calp1 = abs(azi1) ==   90 ? 0 : cos(alp1);
        double phi = lat1 * degree;
        // Ensure cbet1 = +epsilon at poles
        DoublePair bet1 = new DoublePair(_f1 * sin(phi), abs(lat1) == 90 ? Geodesic.tiny : cos(phi));
        Geodesic.sinCosNorm(bet1);
        double sbet1 = bet1.s, cbet1 = bet1.c;
        _dn1 = sqrt(1 + g._ep2 * sq(sbet1));

        // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
        _salp0 = _salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
        // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
        // is slightly better (consider the case salp1 = 0).
        _calp0 = hypot(_calp1, _salp1 * sbet1);
        // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
        // sig = 0 is nearest northward crossing of equator.
        // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
        // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
        // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
        // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
        // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
        // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
        // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
        _somg1 = _salp0 * sbet1;
        _comg1 = sbet1 != 0 || _calp1 != 0 ? cbet1 * _calp1 : 1;
        DoublePair sig1 = new DoublePair(sbet1, _comg1);
        Geodesic.sinCosNorm(sig1);
        _ssig1 = sig1.s; _csig1 = sig1.c;
        // no need to normalize omg1

        _k2 = sq(_calp0) * g._ep2;
        double eps = _k2 / (2 * (1 + sqrt(1 + _k2)) + _k2);
        if ((_caps & CAP_C1) != 0) {
            _A1m1 = Geodesic.A1m1f(eps);
            Geodesic.C1f(eps, _C1a);
            _B11 = Geodesic.sinCosSeries(true, _ssig1, _csig1, _C1a, nC1_);
            double s = sin(_B11), c = cos(_B11);
            // tau1 = sig1 + B11
            _stau1 = _ssig1 * c + _csig1 * s;
            _ctau1 = _csig1 * c - _ssig1 * s;
            // Not necessary because C1pa reverts C1a
            //    _B11 = -SinCosSeries(true, _stau1, _ctau1, _C1pa, nC1p_);
        }

        if ((_caps & CAP_C1p) != 0)
            Geodesic.C1pf(eps, _C1pa);

        if ((_caps & CAP_C2) != 0) {
            _A2m1 = Geodesic.A2m1f(eps);
            Geodesic.C2f(eps, _C2a);
            _B21 = Geodesic.sinCosSeries(true, _ssig1, _csig1, _C2a, nC2_);
        }

        if ((_caps & CAP_C3) != 0) {
            g.C3f(eps, _C3a);
            _A3c = -_f * _salp0 * g.A3f(eps);
            _B31 = Geodesic.sinCosSeries(true, _ssig1, _csig1, _C3a, nC3_ - 1);
        }

        if ((_caps & CAP_C4) != 0) {
            g.C4f(eps, _C4a);
            // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
            _A4 = sq(_a) * _calp0 * _salp0 * g._e2;
            _B41 = Geodesic.sinCosSeries(false, _ssig1, _csig1, _C4a, nC4_);
        }
    }

    public GeodesicLine(Geodesic g, double lat1, double lon1, double azi1) {
        this(g, lat1, lon1, azi1, Geodesic.Mask.ALL);
    }

    public double genPosition(boolean arcmode, double s12_a12, int outmask, Position outPosition) {
        outmask &= _caps & OUT_ALL;
        if (!( init() && (arcmode || ((_caps & DISTANCE_IN & OUT_ALL) != 0)) )) {
            // Uninitialized, or impossible distance calculation requested
            return Double.NaN;
        }

        double sig12, ssig12, csig12, B12 = 0, AB1 = 0;
        if (arcmode) {
            // Interpret s12_a12 as spherical arc length
            sig12 = s12_a12 * degree;
            double s12a = abs(s12_a12);
            s12a -= 180 * floor(s12a / 180);
            ssig12 = s12a ==  0 ? 0 : sin(sig12);
            csig12 = s12a == 90 ? 0 : cos(sig12);
        } else {
            // Interpret s12_a12 as distance
            double tau12 = s12_a12 / (_b * (1 + _A1m1)),
                    s = sin(tau12),
                    c = cos(tau12);
            // tau2 = tau1 + tau12
            B12 = - Geodesic.sinCosSeries(true,
                    _stau1 * c + _ctau1 * s,
                    _ctau1 * c - _stau1 * s,
                    _C1pa, nC1p_);
            sig12 = tau12 - (B12 - _B11);
            ssig12 = sin(sig12); csig12 = cos(sig12);
            if (abs(_f) > 0.01) {
                // Reverted distance series is inaccurate for |f| > 1/100, so correct
                // sig12 with 1 Newton iteration.  The following table shows the
                // approximate maximum error for a = WGS_a() and various f relative to
                // GeodesicExact.
                //     erri = the error in the inverse solution (nm)
                //     errd = the error in the direct solution (series only) (nm)
                //     errda = the error in the direct solution (series + 1 Newton) (nm)
                //
                //       f     erri  errd errda
                //     -1/5    12e6 1.2e9  69e6
                //     -1/10  123e3  12e6 765e3
                //     -1/20   1110 108e3  7155
                //     -1/50  18.63 200.9 27.12
                //     -1/100 18.63 23.78 23.37
                //     -1/150 18.63 21.05 20.26
                //      1/150 22.35 24.73 25.83
                //      1/100 22.35 25.03 25.31
                //      1/50  29.80 231.9 30.44
                //      1/20   5376 146e3  10e3
                //      1/10  829e3  22e6 1.5e6
                //      1/5   157e6 3.8e9 280e6
                double ssig2 = _ssig1 * csig12 + _csig1 * ssig12,
                        csig2 = _csig1 * csig12 - _ssig1 * ssig12;
                B12 = Geodesic.sinCosSeries(true, ssig2, csig2, _C1a, nC1_);
                double serr = (1 + _A1m1) * (sig12 + (B12 - _B11)) - s12_a12 / _b;
                sig12 = sig12 - serr / sqrt(1 + _k2 * sq(ssig2));
                ssig12 = sin(sig12); csig12 = cos(sig12);
                // Update B12 below
            }
        }
        double omg12, lam12, lon12;
        double ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2;
        // sig2 = sig1 + sig12
        ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
        csig2 = _csig1 * csig12 - _ssig1 * ssig12;
        double dn2 = sqrt(1 + _k2 * sq(ssig2));
        if ((outmask & (DISTANCE | REDUCEDLENGTH | GEODESICSCALE)) != 0) {
            if (arcmode || abs(_f) > 0.01)
                B12 = Geodesic.sinCosSeries(true, ssig2, csig2, _C1a, nC1_);
            AB1 = (1 + _A1m1) * (B12 - _B11);
        }
        // sin(bet2) = cos(alp0) * sin(sig2)
        sbet2 = _calp0 * ssig2;
        // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
        cbet2 = hypot(_salp0, _calp0 * csig2);
        if (cbet2 == 0)
            // I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
            cbet2 = csig2 = Geodesic.tiny;
        // tan(omg2) = sin(alp0) * tan(sig2)
        somg2 = _salp0 * ssig2; comg2 = csig2;  // No need to normalize
        // tan(alp0) = cos(sig2)*tan(alp2)
        salp2 = _salp0; calp2 = _calp0 * csig2; // No need to normalize
        // omg12 = omg2 - omg1
        omg12 = atan2(somg2 * _comg1 - comg2 * _somg1,
                comg2 * _comg1 + somg2 * _somg1);
        if ((outmask & DISTANCE) != 0)
            outPosition.s12 = arcmode ? _b * ((1 + _A1m1) * sig12 + AB1) : s12_a12;
        if ((outmask & LONGITUDE) != 0) {
            lam12 = omg12 + _A3c *
                    ( sig12 + (Geodesic.sinCosSeries(true, ssig2, csig2, _C3a, nC3_ - 1)
                            - _B31));
            lon12 = lam12 / degree;
            // Use Math::AngNormalize2 because longitude might have wrapped multiple
            // times.
            lon12 = angNormalize2(lon12);
            outPosition.lon2 = angNormalize(_lon1 + lon12);
        }

        if ((outmask & LATITUDE) != 0)
            outPosition.lat2 = atan2(sbet2, _f1 * cbet2) / degree;

        if ((outmask & AZIMUTH) != 0)
            // minus signs give range [-180, 180). 0- converts -0 to +0.
            outPosition.azi2 = 0 - atan2(-salp2, calp2) / degree;

        if ((outmask & (REDUCEDLENGTH | GEODESICSCALE)) != 0) {
            double B22 = Geodesic.sinCosSeries(true, ssig2, csig2, _C2a, nC2_),
                    AB2 = (1 + _A2m1) * (B22 - _B21),
                    J12 = (_A1m1 - _A2m1) * sig12 + (AB1 - AB2);
            if ((outmask & REDUCEDLENGTH) != 0)
                // Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
                // accurate cancellation in the case of coincident points.
                outPosition.m12 = _b * ((dn2 * (_csig1 * ssig2) - _dn1 * (_ssig1 * csig2))
                        - _csig1 * csig2 * J12);
            if ((outmask & GEODESICSCALE) != 0) {
                double t = _k2 * (ssig2 - _ssig1) * (ssig2 + _ssig1) / (_dn1 + dn2);
                outPosition.M12 = csig12 + (t *  ssig2 -  csig2 * J12) * _ssig1 / _dn1;
                outPosition.M21 = csig12 - (t * _ssig1 - _csig1 * J12) *  ssig2 /  dn2;
            }
        }

        if ((outmask & AREA) != 0) {
            double B42 = Geodesic.sinCosSeries(false, ssig2, csig2, _C4a, nC4_);
            double salp12, calp12;
            if (_calp0 == 0 || _salp0 == 0) {
                // alp12 = alp2 - alp1, used in atan2 so no need to normalized
                salp12 = salp2 * _calp1 - calp2 * _salp1;
                calp12 = calp2 * _calp1 + salp2 * _salp1;
                // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
                // salp12 = -0 and alp12 = -180.  However this depends on the sign being
                // attached to 0 correctly.  The following ensures the correct behavior.
                if (salp12 == 0 && calp12 < 0) {
                    salp12 = Geodesic.tiny * _calp1;
                    calp12 = -1;
                }
            } else {
                // tan(alp) = tan(alp0) * sec(sig)
                // tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
                // = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
                // If csig12 > 0, write
                //   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
                // else
                //   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
                // No need to normalize
                salp12 = _calp0 * _salp0 *
                        (csig12 <= 0 ? _csig1 * (1 - csig12) + ssig12 * _ssig1 :
                                ssig12 * (_csig1 * ssig12 / (1 + csig12) + _ssig1));
                calp12 = sq(_salp0) + sq(_calp0) * _csig1 * csig2;
            }
            outPosition.S12 = _c2 * atan2(salp12, calp12) + _A4 * (B42 - _B41);
        }
        return arcmode ? s12_a12 : sig12 / degree;
    }

    public LatLon position(double s12) {
        Position p = new Position();
        genPosition(false, s12, LATITUDE | LONGITUDE, p);
        return new LatLon(p.lat2, p.lon2);
    }

    public LatLon arcPosition(double a12) {
        Position p = new Position();
        genPosition(true, a12, LATITUDE | LONGITUDE, p);
        return new LatLon(p.lat2, p.lon2);
    }

    public boolean init() {
        return _caps != 0;
    }

    public double latitude() {
        return init() ? _lat1 : Double.NaN;
    }
    public double longitude() {
        return init() ? _lon1 : Double.NaN;
    }
    public double azimuth() {
        return init() ? _azi1 : Double.NaN;
    }
    public double equitorialAzimuth() {
        return init() ? atan2(_salp0, _calp0) / degree : Double.NaN;
    }
    public double equitorialArc() {
        return init() ? atan2(_ssig1, _csig1) / degree : Double.NaN;
    }
    public double majorRadius() {
        return init() ? _a : Double.NaN;
    }
    public double flattening() {
        return init() ? _f : Double.NaN;
    }
    /** @deprecated */
    public double inverseFlattening() {
        return init() ? 1/_f : Double.NaN;
    }
    public int capabilities() {
        return _caps;
    }
    public boolean testCaps(int testcaps) {
        testcaps &= OUT_ALL;
        return (_caps & testcaps) == testcaps;
    }
}
