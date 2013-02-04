/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic;

public final class DoublePair {
    public double s, c;

    public DoublePair(double s, double c) {
        this.s = s;
        this.c = c;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        temp = s != +0.0d ? Double.doubleToLongBits(s) : 0L;
        result = (int)(temp ^ (temp >>> 32));
        temp = c != +0.0d ? Double.doubleToLongBits(c) : 0L;
        result = 31 * result + (int)(temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null || !(obj instanceof DoublePair))
            return false;
        DoublePair o = (DoublePair)obj;
        return s == o.s && c == o.c;
    }
}
