/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic.example;

import com.github.mlc.jgeographic.DMS;
import com.github.mlc.jgeographic.DecodeResult;

public class DMSExample {
    public static void main(String[] args) {
        String dms = "30d14'45.6\"S";
        DecodeResult dr = DMS.decode(dms);
        System.out.println(dr.flag + " " + dr.n);

        double angle = -30.245715;
        System.out.println(DMS.encode(angle, 6, DMS.Flag.LATITUDE));
    }
}
