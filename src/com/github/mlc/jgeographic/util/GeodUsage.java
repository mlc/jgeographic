/*
 * Copyright © Charles Karney (2009-2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * Java port Copyright 2013 Mike Castleman
 * http://github.com/mlc/jgeographic
 */

package com.github.mlc.jgeographic.util;

import java.io.PrintStream;

class GeodUsage {
    static void usage(int retval, boolean brief) {
        PrintStream p = (retval == 0 ? System.out : System.err);
        if (brief) {
            p.print("Usage:\n" +
"    Geod [ -i | -l lat1 lon1 azi1 ] [ -a ] [ -e a f ] [ -d | -: ] [ -b ] [\n" +
"    -f ] [ -p prec ] [ -E ] [ --comment-delimiter commentdelim ] [\n" +
"    --version | -h | --help ] [ --input-file infile | --input-string\n" +
"    instring ] [ --line-separator linesep ] [ --output-file outfile ]\n" +
"\n" +
"For full documentation type:\n" +
"    Geod --help\n" +
"or visit:\n" +
"    http://geographiclib.sf.net/html/Geod.1.html\n");
        } else {
            p.print( "Man page:\n" +
"NAME\n" +
"       Geod -- perform geodesic calculations\n" +
"\n" +
"SYNOPSIS\n" +
"       Geod [ -i | -l lat1 lon1 azi1 ] [ -a ] [ -e a f ] [ -d | -: ] [ -b ] [\n" +
"       -f ] [ -p prec ] [ -E ] [ --comment-delimiter commentdelim ] [\n" +
"       --version | -h | --help ] [ --input-file infile | --input-string\n" +
"       instring ] [ --line-separator linesep ] [ --output-file outfile ]\n" +
"\n" +
"DESCRIPTION\n" +
"       The shortest path between two points on the ellipsoid at (lat1, lon1)\n" +
"       and (lat2, lon2) is called the geodesic.  Its length is s12 and the\n" +
"       geodesic from point 1 to point 2 has forward azimuths azi1 and azi2 at\n" +
"       the two end points.\n" +
"\n" +
"       Geod operates in one of three modes:\n" +
"\n" +
"       1.  By default, Geod accepts lines on the standard input containing\n" +
"           lat1 lon1 azi1 s12 and prints lat2 lon2 azi2 on standard output.\n" +
"           This is the direct geodesic calculation.\n" +
"\n" +
"       2.  Command line arguments -l lat1 lon1 azi1 specify a geodesic line.\n" +
"           Geod then accepts a sequence of s12 values (one per line) on\n" +
"           standard input and prints lat2 lon2 azi2 for each.  This generates\n" +
"           a sequence of points on a single geodesic.\n" +
"\n" +
"       3.  With the -i command line argument, Geod performs the inverse\n" +
"           geodesic calculation.  It reads lines containing lat1 lon1 lat2\n" +
"           lon2 and prints the corresponding values of azi1 azi2 s12.\n" +
"\n" +
"OPTIONS\n" +
"       -i  perform an inverse geodesic calculation (see 3 above).\n" +
"\n" +
"       -l  line mode (see 2 above); generate a sequence of points along the\n" +
"           geodesic specified by lat1 lon1 azi1.\n" +
"\n" +
"       -a  arc mode; on input and output s12 is replaced by a12 the arc length\n" +
"           (in degrees) on the auxiliary sphere.  See \"AUXILIARY SPHERE\".\n" +
"\n" +
"       -e  specify the ellipsoid via a f; the equatorial radius is a and the\n" +
"           flattening is f.  Setting f = 0 results in a sphere.  Specify f < 0\n" +
"           for a prolate ellipsoid.  A simple fraction, e.g., 1/297, is\n" +
"           allowed for f.  (Also, if f > 1, the flattening is set to 1/f.)  By\n" +
"           default, the WGS84 ellipsoid is used, a = 6378137 m, f =\n" +
"           1/298.257223563.\n" +
"\n" +
"       -d  output angles as degrees, minutes, seconds instead of decimal\n" +
"           degrees.\n" +
"\n" +
"       -:  like -d, except use : as a separator instead of the d, ', and \"\n" +
"           delimiters.\n" +
"\n" +
"       -b  report the back azimuth at point 2 instead of the forward azimuth.\n" +
"\n" +
"       -f  full output; each line of output consists of 12 quantities: lat1\n" +
"           lon1 azi1 lat2 lon2 azi2 s12 a12 m12 M12 M21 S12.  a12 is described\n" +
"           in \"AUXILIARY SPHERE\".  The four quantities m12, M12, M21, and S12\n" +
"           are described in \"ADDITIONAL QUANTITIES\".\n" +
"\n" +
"       -p  set the output precision to prec (default 3); prec is the precision\n" +
"           relative to 1 m.  See PRECISION.\n" +
"\n" +
"       -E  use \"exact\" algorithms (based on elliptic integrals) for the\n" +
"           geodesic calculations.  These are more accurate than the (default)\n" +
"           series expansions for |f| > 0.02.\n" +
"\n" +
"       --comment-delimiter\n" +
"           set the comment delimiter to commentdelim (e.g., \"#\" or \"//\").  If\n" +
"           set, the input lines will be scanned for this delimiter and, if\n" +
"           found, the delimiter and the rest of the line will be removed prior\n" +
"           to processing and subsequently appended to the output line\n" +
"           (separated by a space).\n" +
"\n" +
"       --version\n" +
"           print version and exit.\n" +
"\n" +
"       -h  print usage and exit.\n" +
"\n" +
"       --help\n" +
"           print full documentation and exit.\n" +
"\n" +
"       --input-file\n" +
"           read input from the file infile instead of from standard input; a\n" +
"           file name of \"-\" stands for standard input.\n" +
"\n" +
"       --input-string\n" +
"           read input from the string instring instead of from standard input.\n" +
"           All occurrences of the line separator character (default is a\n" +
"           semicolon) in instring are converted to newlines before the reading\n" +
"           begins.\n" +
"\n" +
"       --line-separator\n" +
"           set the line separator character to linesep.  By default this is a\n" +
"           semicolon.\n" +
"\n" +
"       --output-file\n" +
"           write output to the file outfile instead of to standard output; a\n" +
"           file name of \"-\" stands for standard output.\n" +
"\n" +
"INPUT\n" +
"       Geod measures all angles in degrees and all lengths (s12) in meters.\n" +
"       On input angles (latitude, longitude, azimuth, arc length) can be as\n" +
"       decimal degrees or degrees (d), minutes ('), seconds (\").  A decimal\n" +
"       point can only appear in the least significant component and the\n" +
"       designator (d, ', or \") for this component is optional; thus \"40d30\",\n" +
"       \"40d30'\", \"40.5d\", and 40.5 are all equivalent.  By default, latitude\n" +
"       precedes longitude for each point; however on input either may be given\n" +
"       first by appending (or prepending) N or S to the latitude and E or W to\n" +
"       the longitude.  Azimuths are measured clockwise from north; however\n" +
"       this may be overridden with E or W.\n" +
"\n" +
"       See the \"QUOTING\" section of GeoConvert(1) for how to quote the DMS\n" +
"       designators ' and \".\n" +
"\n" +
"AUXILIARY SPHERE\n" +
"       Geodesics on the ellipsoid can be transferred to the auxiliary sphere\n" +
"       on which the distance is measured in terms of the arc length a12\n" +
"       (measured in degrees) instead of s12.  In terms of a12, 180 degrees is\n" +
"       the distance from one equator crossing to the next or from the minimum\n" +
"       latitude to the maximum latitude.  Geodesics with a12 > 180 degrees do\n" +
"       not correspond to shortest paths.  With the -a flag, s12 (on both input\n" +
"       and output) is replaced by a12.  The -a flag does not affect the full\n" +
"       output given by the -f flag (which always includes both s12 and a12).\n" +
"\n" +
"ADDITIONAL QUANTITIES\n" +
"       The -f flag reports four additional quantities.\n" +
"\n" +
"       The reduced length of the geodesic, m12, is defined such that if the\n" +
"       initial azimuth is perturbed by dazi1 (radians) then the second point\n" +
"       is displaced by m12 dazi1 in the direction perpendicular to the\n" +
"       geodesic.  m12 is given in meters.  On a curved surface the reduced\n" +
"       length obeys a symmetry relation, m12 + m21 = 0.  On a flat surface, we\n" +
"       have m12 = s12.\n" +
"\n" +
"       M12 and M21 are geodesic scales.  If two geodesics are parallel at\n" +
"       point 1 and separated by a small distance dt, then they are separated\n" +
"       by a distance M12 dt at point 2.  M21 is defined similarly (with the\n" +
"       geodesics being parallel to one another at point 2).  M12 and M21 are\n" +
"       dimensionless quantities.  On a flat surface, we have M12 = M21 = 1.\n" +
"\n" +
"       If points 1, 2, and 3 lie on a single geodesic, then the following\n" +
"       addition rules hold, s13 = s12 + s23, a13 = a12 + a23, S13 = S12 + S23,\n" +
"       m13 = m12 M23 + m23 M21, M13 = M12 M23 - (1 - M12 M21) m23 / m12, and\n" +
"       M31 = M32 M21 - (1 - M23 M32) m12 / m23.\n" +
"\n" +
"       Finally, S12 is the area between the geodesic from point 1 to point 2\n" +
"       and the equator; i.e., it is the area, measured counter-clockwise, of\n" +
"       the quadrilateral with corners (lat1,lon1), (0,lon1), (0,lon2), and\n" +
"       (lat2,lon2).  It is given in meters^2.\n" +
"\n" +
"PRECISION\n" +
"       prec gives precision of the output with prec = 0 giving 1 m precision,\n" +
"       prec = 3 giving 1 mm precision, etc.  prec is the number of digits\n" +
"       after the decimal point for lengths.  For decimal degrees, the number\n" +
"       of digits after the decimal point is 5 + prec.  For DMS (degree,\n" +
"       minute, seconds) output, the number of digits after the decimal point\n" +
"       in the seconds component is 1 + prec.  The minimum value of prec is 0\n" +
"       and the maximum is 10.\n" +
"\n" +
"ERRORS\n" +
"       An illegal line of input will print an error message to standard output\n" +
"       beginning with \"ERROR:\" and causes Geod to return an exit code of 1.\n" +
"       However, an error does not cause Geod to terminate; following lines\n" +
"       will be converted.\n" +
"\n" +
"ACCURACY\n" +
"       Using the (default) series solution, Geod is accurate to about 15 nm\n" +
"       (15 nanometers) for the WGS84 ellipsoid.  The approximate maximum error\n" +
"       (expressed as a distance) for an ellipsoid with the same major radius\n" +
"       as the WGS84 ellipsoid and different values of the flattening is\n" +
"\n" +
"          |f|     error\n" +
"          0.01    25 nm\n" +
"          0.02    30 nm\n" +
"          0.05    10 um\n" +
"          0.1    1.5 mm\n" +
"          0.2    300 mm\n" +
"\n" +
"       If -E is specified, Geod is accurate to about 40 nm (40 nanometers) for\n" +
"       the WGS84 ellipsoid.  The approximate maximum error (expressed as a\n" +
"       distance) for an ellipsoid with a quarter meridian of 10000 km and\n" +
"       different values of the a/b = 1 - f is\n" +
"\n" +
"          1-f    error (nm)\n" +
"          1/128   387\n" +
"          1/64    345\n" +
"          1/32    269\n" +
"          1/16    210\n" +
"          1/8     115\n" +
"          1/4      69\n" +
"          1/2      36\n" +
"            1      15\n" +
"            2      25\n" +
"            4      96\n" +
"            8     318\n" +
"           16     985\n" +
"           32    2352\n" +
"           64    6008\n" +
"          128   19024\n" +
"\n" +
"MULTIPLE SOLUTIONS\n" +
"       The shortest distance returned for the inverse problem is (obviously)\n" +
"       uniquely defined.  However, in a few special cases there are multiple\n" +
"       azimuths which yield the same shortest distance.  Here is a catalog of\n" +
"       those cases:\n" +
"\n" +
"       lat1 = -lat2 (with neither at a pole)\n" +
"           If azi1 = azi2, the geodesic is unique.  Otherwise there are two\n" +
"           geodesics and the second one is obtained by setting [azi1,azi2] =\n" +
"           [azi2,azi1], [M12,M21] = [M21,M12], S12 = -S12.  (This occurs when\n" +
"           the longitude difference is near +/-180 for oblate ellipsoids.)\n" +
"\n" +
"       lon2 = lon1 +/- 180 (with neither at a pole)\n" +
"           If azi1 = 0 or +/-180, the geodesic is unique.  Otherwise there are\n" +
"           two geodesics and the second one is obtained by setting [azi1,azi2]\n" +
"           = [-azi1,-azi2], S12 = -S12.  (This occurs when the lat2 is near\n" +
"           -lat1 for prolate ellipsoids.)\n" +
"\n" +
"       Points 1 and 2 at opposite poles\n" +
"           There are infinitely many geodesics which can be generated by\n" +
"           setting [azi1,azi2] = [azi1,azi2] + [d,-d], for arbitrary d.  (For\n" +
"           spheres, this prescription applies when points 1 and 2 are\n" +
"           antipodal.)\n" +
"\n" +
"       s12 = 0 (coincident points)\n" +
"           There are infinitely many geodesics which can be generated by\n" +
"           setting [azi1,azi2] = [azi1,azi2] + [d,d], for arbitrary d.\n" +
"\n" +
"EXAMPLES\n" +
"       Route from JFK Airport to Singapore Changi Airport:\n" +
"\n" +
"          echo 40:38:23N 073:46:44W 01:21:33N 103:59:22E |\n" +
"          Geod -i -: -p 0\n" +
"\n" +
"          003:18:29.9 177:29:09.2 15347628\n" +
"\n" +
"       Waypoints on the route at intervals of 2000km:\n" +
"\n" +
"          for ((i = 0; i <= 16; i += 2)); do echo ${i}000000;done |\n" +
"          Geod -l 40:38:23N 073:46:44W 003:18:29.9 -: -p 0\n" +
"\n" +
"          40:38:23.0N 073:46:44.0W 003:18:29.9\n" +
"          58:34:45.1N 071:49:36.7W 004:48:48.8\n" +
"          76:22:28.4N 065:32:17.8W 010:41:38.4\n" +
"          84:50:28.0N 075:04:39.2E 150:55:00.9\n" +
"          67:26:20.3N 098:00:51.2E 173:27:20.3\n" +
"          49:33:03.2N 101:06:52.6E 176:07:54.3\n" +
"          31:34:16.5N 102:30:46.3E 177:03:08.4\n" +
"          13:31:56.0N 103:26:50.7E 177:24:55.0\n" +
"          04:32:05.7S 104:14:48.7E 177:28:43.6\n" +
"\n" +
"WARNING\n" +
"       The proj library has a utility geod which offers similar functionality\n" +
"       (but which uses less accurate algorithms).  Geod and geod are distinct\n" +
"       names on most operating systems, but can be confused on Windows\n" +
"       systems.  If you have proj installed on your Windows system, then in\n" +
"       order to use Geod you will need to use the full pathname for the\n" +
"       executable or ensure that your \"PATH\" includes the GeographicLib binary\n" +
"       directory before the proj directory.\n" +
"\n" +
"SEE ALSO\n" +
"       GeoConvert(1).  The algorithms are described in C. F. F. Karney,\n" +
"       Algorithms for geodesics, J. Geodesy 87, 43-55 (2013); DOI:\n" +
"       http://dx.doi.org/10.1007/s00190-012-0578-z\n" +
"       <http://dx.doi.org/10.1007/s00190-012-0578-z>; addenda:\n" +
"       http://geographiclib.sf.net/geod-addenda.html\n" +
"       <http://geographiclib.sf.net/geod-addenda.html>.\n" +
"\n" +
"AUTHOR\n" +
"       Geod was written by Charles Karney.\n" +
"\n" +
"HISTORY\n" +
"       Geod was added to GeographicLib, <http://geographiclib.sf.net>, in\n" +
"       2009-03.\n");
        }
        System.exit(retval);
    }
}