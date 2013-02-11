/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic.util;

import java.io.*;
import java.lang.Math;
import com.github.mlc.jgeographic.*;

public class Geod {
    private static String latLonString(double lat, double lon, int prec, boolean dms, char dmssep) {
          return dms ?
                  DMS.encode(lat, prec + 5, DMS.Flag.LATITUDE, dmssep) + " " +
                          DMS.encode(lon, prec + 5, DMS.Flag.LONGITUDE, dmssep) :
                  DMS.encode(lat, prec + 5, DMS.Flag.NUMBER) + " " +
                          DMS.encode(lon, prec + 5, DMS.Flag.NUMBER);
    }

    private static String azimuthString(double azi, int prec, boolean  dms, char dmssep) {
        return dms ? DMS.encode(azi, prec + 5, DMS.Flag.AZIMUTH, dmssep) :
                DMS.encode(azi >= 180 ? azi - 360 : azi, prec + 5, DMS.Flag.NUMBER);
    }

    private static String distanceStrings(double s12, double a12,
                                          boolean full, boolean arcmode, int prec, boolean dms) {
        StringBuilder s = new StringBuilder();
        if (full || !arcmode) {
            s.append(Utility.str(s12, prec));
        }
        if (full)
            s.append(" ");
        if (full || arcmode)
            s.append(DMS.encode(a12, prec + 5, dms ? DMS.Flag.NONE : DMS.Flag.NUMBER));

        return s.toString();
    }

    private static double readDistance(String s, boolean arcmode) {
        return arcmode ? DMS.decodeAngle(s) : Utility.num(s);
    }

    public static void main(String[] args) throws IOException {
        boolean linecalc = false, inverse = false, arcmode = false,
                dms = false, full = false, exact = false;
        double a = Constants.WGS84_a,
                f = Constants.WGS84_f;
        double lat1 = 0.0, lon1 = 0.0, azi1 = 0.0, a12, s12;
        double azi2sense = 0;
        Position position = new Position();
        InversePosition inversePosition = new InversePosition();
        int prec = 3;
        String istring = null, ifile = null, ofile = null, cdelim = null;
        char lsep=';', dmssep='\0';

        for (int m = 0, argc = args.length; m < argc; ++m) {
            String arg = args[m];
            if (arg.equals("-i")) {
                inverse = true;
                linecalc = false;
            } else if (arg.equals("-a")) {
                arcmode = true;
            } else if (arg.equals("-l")) {
                inverse = false;
                linecalc = true;
                if (m + 3 >= argc) {
                    GeodUsage.usage(1, true);
                }
                try {
                    LatLon latLon = DMS.decodeLatLon(args[m + 1], args[m + 2]);
                    lat1 = latLon.lat;
                    lon1 = latLon.lon;
                    azi1 = DMS.decodeAngle(args[m+3]);
                } catch (Exception ex) {
                    System.err.println("Error decoding arguments of -l: " + ex.toString());
                    System.exit(1);
                }
                m += 3;
            } else if (arg.equals("-e")) {
                if (m + 2 >= argc) {
                    GeodUsage.usage(1, true);
                }
                try {
                    a = Utility.num(args[m + 1]);
                    f = Utility.fract(args[m+2]);
                } catch (Exception ex) {
                    System.err.println("Error decoding arguments of -e: " + ex.toString());
                    System.exit(1);
                }
                m += 2;
            } else if (arg.equals("-d")) {
                dms = true;
                dmssep = '\0';
            } else if (arg.equals("-:")) {
                dms = true;
                dmssep = ':';
            } else if (arg.equals("-b"))
                azi2sense = 180;
            else if (arg.equals("-f"))
                full = true;
            else if (arg.equals("-p")) {
                if (++m == argc)
                    GeodUsage.usage(1, true);
                try {
                    prec = Integer.parseInt(args[m], 10);
                } catch (Exception ex) {
                    System.err.println("Precision " + args[m] + " is not a number");
                    System.exit(1);
                }
            } else if (arg.equals("-E"))
                exact = true;
            else if (arg.equals("--input-string")) {
                if (++m == argc) GeodUsage.usage(1, true);
                istring = args[m];
            } else if (arg.equals("--input-file")) {
                if (++m == argc) GeodUsage.usage(1, true);
                ifile = args[m];
            } else if (arg.equals("--output-file")) {
                if (++m == argc) GeodUsage.usage(1, true);
                ofile = args[m];
            } else if (arg.equals("--line-separator")) {
                if (++m == argc) GeodUsage.usage(1, true);
                if (args[m].length() != 1) {
                    System.err.println("Line separator must be a single character");
                    System.exit(1);
                }
                lsep = args[m].charAt(0);
            } else if (arg.equals("--version")) {
                System.out.println("Geod from jgeodsic");
                return;
            } else
                GeodUsage.usage(!(arg.equals("-h") || arg.equals("--help")) ? 1 : 0, !arg.equals("--help"));
        }
        if (ifile != null && istring != null) {
            System.err.println("Cannot specify --input-string and --input-file together");
            System.exit(1);
        }
        if ("-".equals(ifile))
            ifile = null;
        BufferedReader input;
        try {
            if (ifile != null) {
                input = new BufferedReader(new FileReader(ifile));
            } else if (istring != null) {
                input = new BufferedReader(new StringReader(istring.replaceAll(String.valueOf(lsep), "\n")));
            } else {
                input = new BufferedReader(new InputStreamReader(System.in));
            }
        } catch (IOException ex) {
            System.err.println("Can not open input: " + ex.toString());
            System.exit(1);
            return; // shut up java
        }

        PrintWriter output;
        try {
            if (ofile != null) {
                output = new PrintWriter(new FileWriter(ofile));
            } else {
                output = new PrintWriter(System.out);
            }
        } catch (IOException ex) {
            System.err.println("Can not open output: " + ex.toString());
            System.exit(1);
            return; // shut up java
        }

        Geodesic geod = new Geodesic(a, f);
        GeodesicLine l;
        if (linecalc) {
            l = geod.line(lat1, lon1, azi1);
        } else {
            l = null; // shut up java
        }
        prec = Math.min(10, Math.max(0, prec));
        String s;
        int retval = 0;
        while ((s = input.readLine()) != null) {
            try {
                String eol = "";
                if (cdelim != null) {
                    int m = s.indexOf(cdelim);
                    if (m >= 0) {
                        eol = " " + s.substring(m);
                        s = s.substring(0, m);
                    }
                }
                String[] components = s.split("\\s+");
                if (inverse) {
                    if (components.length != 4) {
                        throw new GeographicException("invalid input " + s);
                    }
                    LatLon ll1 = DMS.decodeLatLon(components[0], components[1]),
                            ll2 = DMS.decodeLatLon(components[2], components[3]);
                    a12 = geod.genInverse(ll1.lat, ll1.lon, ll2.lat, ll2.lon, Geodesic.Mask.ALL, inversePosition);
                    if (full) {
                        output.print(latLonString(ll1.lat, ll1.lon, prec, dms, dmssep));
                        output.print(' ');
                    }
                    output.print(azimuthString(inversePosition.azi1, prec, dms, dmssep));
                    output.print(' ');
                    if (full) {
                        output.print(latLonString(ll2.lat, ll2.lon, prec, dms, dmssep));
                        output.print(' ');
                    }
                    output.print(azimuthString(inversePosition.azi2 + azi2sense, prec, dms, dmssep));
                    output.print(' ');
                    output.print(distanceStrings(inversePosition.s12, a12, full, arcmode, prec, dms));
                    if (full) {
                        output.print(' '); output.print(Utility.str(inversePosition.m12, prec));
                        output.print(' '); output.print(Utility.str(inversePosition.M12, prec + 7));
                        output.print(' '); output.print(Utility.str(inversePosition.M21, prec + 7));
                        output.print(' '); output.print(Utility.str(inversePosition.S12, Math.max(prec - 7 , 0)));
                    }
                    output.println(eol);
                } else {
                    if (linecalc) {
                        if (components.length != 1) {
                            throw new GeographicException("invalid input " + s);
                        }
                        s12 = readDistance(components[0], arcmode);
                        a12 = l.genPosition(arcmode, s12, Geodesic.Mask.ALL, position);
                    } else {
                        if (components.length != 4) {
                            throw new GeographicException("invalid input " + s);
                        }
                        LatLon latLon = DMS.decodeLatLon(components[0], components[1]);
                        lat1 = latLon.lat; lon1 = latLon.lon;
                        azi1 = DMS.decodeAzimuth(components[2]);
                        s12 = readDistance(components[3], arcmode);
                        a12 = geod.genDirect(lat1, lon1, azi1, arcmode, s12, Geodesic.Mask.ALL, position);
                    }
                    if (arcmode) {
                        double tmp = s12; s12 = a12; a12 = tmp;
                    }
                    if (full) {
                        output.print(latLonString(lat1, lon1, prec, dms, dmssep));
                        output.print(' ');
                        output.print(azimuthString(azi1, prec, dms, dmssep));
                        output.print(' ');
                    }
                    output.print(latLonString(position.lat2, position.lon2, prec, dms, dmssep));
                    output.print(' ');
                    output.print(azimuthString(position.azi2 + azi2sense, prec, dms, dmssep));
                    if (full) {
                        output.print(' ');
                        output.print(distanceStrings(s12, a12, full, arcmode, prec, dms));
                        output.print(' '); output.print(Utility.str(position.m12, prec));
                        output.print(' '); output.print(Utility.str(position.M12, prec + 7));
                        output.print(' '); output.print(Utility.str(position.M21, prec + 7));
                        output.print(' '); output.print(Utility.str(position.S12, Math.max(prec - 7, 0)));
                    }
                    output.println(eol);
                }
            } catch (Exception ex) {
                output.println("ERROR: " + ex.toString());
                retval = 1;
            }
        }
        output.close();
        System.exit(retval);
    }
}
