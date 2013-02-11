/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic.util;

import java.text.NumberFormat;
import java.util.Locale;
import com.github.mlc.jgeographic.GeographicException;

public class Utility {
    private Utility() {}
    public static String str(double x, int p) {
        if (Double.isInfinite(x)) {
            return (x < 0 ? "-inf" : "inf");
        } else if (Double.isNaN(x)) {
            return "nan";
        }
        if (p >= 0) {
            NumberFormat nf = NumberFormat.getInstance(Locale.US);
            nf.setGroupingUsed(false);
            nf.setMaximumFractionDigits(p);
            nf.setMinimumFractionDigits(p);
            return nf.format(x);
        } else {
            return Double.toString(x);
        }
    }

    public static double num(String s) {
        String errmsg;
        try {
            NumberFormat nf = NumberFormat.getInstance(Locale.US);
            Number x = nf.parse(s);
            return x.doubleValue();
        } catch (Exception e) {
            errmsg = e.getMessage();
        }

        double x = nummatch(s);
        if (x == 0)
            throw new GeographicException(errmsg);
        return x;
    }

    public static double nummatch(String s) {
        if (s.length() < 3)
            return 0;
        String t = s.toUpperCase(Locale.US);
        int sign = t.charAt(0) == '-' ? -1 : 1;
        int p0 = t.charAt(0) == '-' || t.charAt(0) == '+' ? 1 : 0;
        int p1 = find_last_not_of(t, '0');
        if (p1 == -1 || p1 + 1 < p0 + 3)
            return 0;
        // Strip off sign and trailing 0s
        t = t.substring(p0, p1 + 1);
        if (t.equals("NAN") || t.equals("1.#QNAN") || t.equals("1.#SNAN") || t.equals("1.#IND") ||
                t.equals("1.#R"))
            return Double.NaN;
        else if (t.equals("INF") || t.equals("1.#INF"))
            return sign < 1 ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        return 0;
    }

    private static int find_last_not_of(String haystack, char needle) {
        for (int pos = haystack.length() - 1; pos >= 0; --pos) {
            if (haystack.charAt(pos) != needle)
                return pos;
        }
        return -1;
    }

    public static double fract(String s) {
        int delim = s.indexOf('/');
        return !(delim >= 1 && delim + 2 <= s.length()) ?
                num(s) :
                // delim in [1, size() - 2]
                num(s.substring(0, delim)) / num(s.substring(delim + 1));
    }
}