/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic;

import java.text.NumberFormat;
import java.util.Locale;

import static com.github.mlc.jgeographic.Math.angNormalize;
import static java.lang.Math.*;

public class DMS {
    private DMS() { }

    private static final char[] hemispheres = {'S', 'N', 'W', 'E'};
    private static final char[] signs = {'-', '+'};
    private static final char[] digits = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
    private static final char[] dmsindicators = {'D', '\'', '"', ':'};
    private static final String[] components = { "degrees", "minutes", "seconds" };

    public static enum Flag {
        NONE, LATITUDE, LONGITUDE, AZIMUTH, NUMBER
    }

    public static enum Component {
        DEGREE, MINUTE, SECOND
    }

    public static double decode(double d, double m, double s) {
        return d + (m + s/60.0)/60.0;
    }
    public static double decode(double d, double m) {
        return decode(d, m, 0.0);
    }
    public static double decode(double d) {
        return decode(d, 0.0, 0.0);
    }

    public static DecodeResult decode(String dms) {
        String errormsg = null;
        char[] dmsa = dms.toCharArray();
        int end = dmsa.length;
        for (int i = 0; i < end; ++i) {
            switch(dmsa[i]) {
            case '\u00B0': // degree symbol
            case '\u00BA': // alt symbol
            case '\u2070': // sup zero
            case '\u02DA': // ring above
                dmsa[i] = 'd';
                break;
            case '\u2032': // prime
            case '\u00B4': // acute accent
            case '\u2019': // right single quote
                dmsa[i] = '\'';
                break;
            case '\u2033': // double prime
            case '\u201D': // right double quote
                dmsa[i] = '"';
                break;
            }
            if (i > 0 && dmsa[i] == '\'' && dmsa[i-1] == '\'') {
                i--;
                dmsa[i] = '"';
                System.arraycopy(dmsa, i+2, dmsa, i+1, (end - i - 1));
                end--;
            }
        }
        do {
            int sign = 1;
            int beg = 0;
            while(beg < end && Character.isWhitespace(dmsa[beg]))
                ++beg;
            while (beg < end && Character.isWhitespace(dmsa[end-1]))
                --end;
            Flag ind1 = Flag.NONE;
            int k;
            if (end > beg && (k = lookup(hemispheres, dmsa[beg])) >= 0) {
                ind1 = (k / 2 != 0) ? Flag.LONGITUDE : Flag.LATITUDE;
                sign = (k % 2 != 0) ? 1 : -1;
                ++beg;
            }
            if (end > beg && (k = lookup(hemispheres, dmsa[end-1])) >= 0) {
                if (k >= 0) {
                    if (ind1 != Flag.NONE) {
                        if (Character.toUpperCase(dmsa[beg - 1]) == Character.toUpperCase(dmsa[end - 1]))
                            errormsg = "Repeated hemisphere indicators "
                                    + dmsa[beg - 1]
                                    + " in " + new String(dmsa, beg - 1, end - beg + 1);
                        else
                            errormsg = "Contradictory hemisphere indicators "
                                    + dmsa[beg - 1] + " and "
                                    + dmsa[end - 1] + " in "
                                    + new String(dmsa, beg - 1, end - beg + 1);
                        break;
                    }
                    ind1 = (k / 2 != 0) ? Flag.LONGITUDE : Flag.LATITUDE;
                    sign = (k % 2 != 0) ? 1 : -1;
                    --end;
                }
            }
            if (end > beg && (k = lookup(signs, dmsa[beg])) >= 0) {
                if (k >= 0) {
                    sign *= (k != 0) ? 1 : -1;
                    ++beg;
                }
            }
            if (end == beg) {
                errormsg = "Empty or incomplete DMS string " + new String(dmsa);
            }
            double[] ipieces = new double[3];
            double[] fpieces = new double[3];
            int npiece = 0;
            double icurrent = 0;
            double fcurrent = 0;
            int ncurrent = 0, p = beg;
            boolean pointseen = false;
            int digcount = 0, intcount = 0;
            while (p < end) {
                char x = dmsa[p++];
                if ((k = lookup(digits, x)) >= 0) {
                    ++ncurrent;
                    if (digcount > 0)
                        ++digcount;         // Count of decimal digits
                    else {
                        icurrent = 10 * icurrent + k;
                        ++intcount;
                    }
                } else if (x == '.') {
                    if (pointseen) {
                        errormsg = "Multiple decimal points in "
                                + new String(dmsa, beg, end - beg);
                        break;
                    }
                    pointseen = true;
                    digcount = 1;
                } else if ((k = lookup(dmsindicators, x)) >= 0) {
                    if (k >= 3) {
                        if (p == end) {
                            errormsg = "Illegal for : to appear at the end of " +
                                    new String(dmsa, beg, end - beg);
                            break;
                        }
                        k = npiece;
                    }
                    if (k == npiece - 1) {
                        errormsg = "Repeated " + components[k] +
                                " component in " + new String(dmsa, beg, end - beg);
                        break;
                    } else if (k < npiece) {
                        errormsg = components[k] + " component follows "
                                + components[npiece - 1] + " component in "
                                + new String(dmsa, beg, end - beg);
                        break;
                    }
                    if (ncurrent == 0) {
                        errormsg = "Missing numbers in " + components[k] +
                                " component of " + new String(dmsa, beg, end - beg);
                        break;
                    }
                    if (digcount > 1) {
                        String s = new String(dmsa, p - intcount - digcount - 1,
                                intcount + digcount);
                        fcurrent = Double.parseDouble(s);
                        icurrent = 0;
                    }
                    ipieces[k] = icurrent;
                    fpieces[k] = icurrent + fcurrent;
                    if (p < end) {
                        npiece = k + 1;
                        icurrent = fcurrent = 0;
                        ncurrent = digcount = intcount = 0;
                    }
                } else if (lookup(signs, x) >= 0) {
                    errormsg = "Internal sign in DMS string "
                            + new String(dmsa, beg, end - beg);
                    break;
                } else {
                    errormsg = "Illegal character " + x + " in DMS string "
                            + new String(dmsa, beg, end - beg);
                    break;
                }
            }
            if (errormsg != null)
                break;
            if (lookup(dmsindicators, dmsa[p - 1]) < 0) {
                if (npiece >= 3) {
                    errormsg = "Extra text following seconds in DMS string "
                            + new String(dmsa, beg, end - beg);
                    break;
                }
                if (ncurrent == 0) {
                    errormsg = "Missing numbers in trailing component of "
                            + new String(dmsa, beg, end - beg);
                    break;
                }
                if (digcount > 1) {
                    String s = new String(dmsa, p - intcount - digcount,
                            intcount + digcount);
                    fcurrent = Double.parseDouble(s);
                    icurrent = 0;
                }
                ipieces[npiece] = icurrent;
                fpieces[npiece] = icurrent + fcurrent;
            }
            if (pointseen && digcount == 0) {
                errormsg = "Decimal point in non-terminal component of "
                        + new String(dmsa, beg, end - beg);
                break;
            }
            // Note that we accept 59.999999... even though it rounds to 60.
            if (ipieces[1] >= 60) {
                errormsg = "Minutes " + fpieces[1]
                        + " not in range [0, 60)";
                break;
            }
            if (ipieces[2] >= 60) {
                errormsg = "Seconds " + fpieces[2]
                        + " not in range [0, 60)";
                break;
            }
            //ind = ind1;
            // Assume check on range of result is made by calling routine (which
            // might be able to offer a better diagnostic).
            return new DecodeResult(sign * (fpieces[0] + (fpieces[1] + fpieces[2] / 60) / 60), ind1);

        } while(false);
        double val = 0.0;
        try {
            Double.parseDouble(String.valueOf(dmsa));
        } catch (Exception ignore) {
        }
        if (val == 0)
            throw new GeographicException(errormsg);
        return new DecodeResult(val, Flag.NONE);
    }

    public static LatLon decodeLatLon(String stra, String strb, boolean swaplatlong) {
        DecodeResult ra = decode(stra);
        DecodeResult rb = decode(strb);
        double a = ra.n, b = rb.n;
        Flag ia = ra.flag, ib = rb.flag;

        if (ia == Flag.NONE && ib == Flag.NONE) {
            ia = swaplatlong ? Flag.LONGITUDE : Flag.LATITUDE;
            ib = swaplatlong ? Flag.LATITUDE : Flag.LONGITUDE;
        } else if (ia == Flag.NONE) {
            ia = (ib == Flag.LATITUDE) ? Flag.LONGITUDE : Flag.LATITUDE;
        } else if (ib == Flag.NONE) {
            ib = (ia == Flag.LATITUDE) ? Flag.LONGITUDE : Flag.LATITUDE;
        }
        if (ia == ib) {
            throw new GeographicException("Both " + stra + " and "
                    + strb + " interpreted as "
                    + (ia == Flag.LATITUDE ? "latitudes" : "longitudes"));
        }
            double lat1 = ia == Flag.LATITUDE ? a : b,
                    lon1 = ia == Flag.LATITUDE ? b : a;
        if (abs(lat1) > 90)
            throw new GeographicException("Latitude " + lat1
                + "d not in [-90d, 90d]");
        if (lon1 < -540 || lon1 >= 540)
            throw new GeographicException("Longitude " + lon1
                + "d not in [-540d, 540d)");
        lon1 = angNormalize(lon1);
        return new LatLon(lat1, lon1);
    }

    public static LatLon decodeLatLon(String stra, String strb) {
        return decodeLatLon(stra, strb, false);
    }

    public static double decodeAngle(String angstr) {
        DecodeResult r = decode(angstr);
        if (r.flag != Flag.NONE) {
            throw new GeographicException("Arc angle " + angstr
                    + " includes a hemisphere, N/E/W/S");
        }
        return r.n;
    }

    public static double decodeAzimuth(String azistr) {
        DecodeResult r = decode(azistr);
        if (r.flag == Flag.LATITUDE) {
            throw new GeographicException("Azimuth " + azistr
                          + " has a latitude hemisphere, N/S");
        }
        if (r.n < -540 || r.n >= 540)
            throw new GeographicException("Azimuth " + azistr + " not in range [-540d, 540d)");
        return angNormalize(r.n);
    }

    public static String encode(double angle, Component trailing, int prec, Flag ind, char dmssep) {
        if (Double.isInfinite(angle)) {
            return angle < 0 ? "-inf" : "inf";
        } else if (Double.isNaN(angle)) {
            return "nan";
        }
        // 15 - 2 * trailing = ceiling(log10(2^53/90/60^trailing)).
        // This suffices to give full real precision for numbers in [-90,90]
        prec = min(15 - 2 * trailing.ordinal(), prec);
        double scale = 1;
        for (int i = 0; i < trailing.ordinal(); ++i)
            scale *= 60;
        for (int i = 0; i < prec; ++i) {
            scale *= 10;
        }
        if (ind == Flag.AZIMUTH)
            angle -= floor(angle/360) * 360;
        int sign = angle < 0 ? -1 : 1;
        angle *= sign;
        double idegree = floor(angle),
                fdegree = floor((angle - idegree) * scale + 0.5) / scale;
        if (fdegree >= 1) {
            idegree += 1;
            fdegree -= 1;
        }
        double[] pieces = {fdegree, 0, 0};
        for (int i = 1; i <= trailing.ordinal(); ++i) {
            double ip = floor(pieces[i - 1]),
                    fp = pieces[i - 1] - ip;
            pieces[i] = fp * 60;
            pieces[i - 1] = ip;
        }
        pieces[0] += idegree;
        StringBuilder s = new StringBuilder();
        NumberFormat nf = NumberFormat.getInstance(Locale.US);
        if (ind == Flag.NONE && sign < 0)
            s.append('-');
        switch (trailing) {
        case DEGREE:
            if (ind != Flag.NONE) {
                nf.setMinimumIntegerDigits(1 + min(ind.ordinal(), 2));
            }
            nf.setMinimumFractionDigits(prec);
            nf.setMaximumFractionDigits(prec);
            s.append(nf.format(pieces[0])).append(dmssep != '\0' ? dmssep : dmsindicators[0]);
            break;
        default:
            if (ind != Flag.NONE) {
                nf.setMinimumIntegerDigits(1 + min(ind.ordinal(), 2));
            }
            nf.setMaximumFractionDigits(0);
            s.append(nf.format(pieces[0])).append(dmssep != '\0' ? dmssep : Character.toLowerCase(dmsindicators[0]));
            switch (trailing) {
            case MINUTE:
                nf.setMinimumIntegerDigits(2);
                nf.setMinimumFractionDigits(prec);
                nf.setMaximumFractionDigits(prec);
                s.append(nf.format(pieces[1]));
                if (dmssep != '\0')
                    s.append(Character.toLowerCase(dmsindicators[1]));
                break;
            case SECOND:
                nf.setMinimumIntegerDigits(2);
                s.append(nf.format(pieces[1])).append(dmssep != '\0' ? dmssep : Character.toLowerCase(dmsindicators[1]));
                nf.setMinimumFractionDigits(prec);
                nf.setMaximumFractionDigits(prec);
                s.append(nf.format(pieces[2]));
                if (dmssep == '\0')
                    s.append(Character.toLowerCase(dmsindicators[2]));
                break;
            }
        }
        if (ind != Flag.NONE && ind != Flag.AZIMUTH)
            s.append(hemispheres[ind == Flag.LATITUDE ? 0 : 2 + (sign < 0 ? 0 : 1)]);
        return s.toString();
    }

    public static String encode(double angle, Component trailing, int prec, Flag ind) {
        return encode(angle, trailing, prec, ind, '\0');
    }
    public static String encode(double angle, Component trailing, int prec) {
        return encode(angle, trailing, prec, Flag.NONE, '\0');
    }
    public static String encode(double angle, int prec, Flag ind, char dmssep) {
        return encode(angle,
                prec < 2 ? Component.DEGREE : (prec < 4 ? Component.MINUTE : Component.SECOND),
                prec < 2 ? prec : (prec < 4 ? prec - 2 : prec - 4),
                ind, dmssep);
    }
    public static String encode(double angle, int prec, Flag ind) {
        return encode(angle, prec, ind, '\0');
    }
    public static String encode(double angle, int prec) {
        return encode(angle, prec, Flag.NONE, '\0');
    }
    private static int lookup(char[] haystack, char needle) {
        needle = Character.toUpperCase(needle);
        for (int i = 0, len = haystack.length; i < len; ++i) {
            if (haystack[i] == needle)
                return i;
        }
        return -1;
    }
}
