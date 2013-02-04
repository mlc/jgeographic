/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic;

public class DecodeResult {
    public final double n;
    public final DMS.Flag flag;

    public DecodeResult(double n, DMS.Flag flag) {
        this.n = n;
        this.flag = flag;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof DecodeResult)) return false;

        DecodeResult that = (DecodeResult)o;

        if (Double.compare(that.n, n) != 0) return false;
        if (flag != that.flag) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        temp = n != +0.0d ? Double.doubleToLongBits(n) : 0L;
        result = (int)(temp ^ (temp >>> 32));
        result = 31 * result + flag.hashCode();
        return result;
    }
}
