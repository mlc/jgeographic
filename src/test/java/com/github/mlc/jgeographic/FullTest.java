package com.github.mlc.jgeographic;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.zip.GZIPInputStream;
import com.google.common.collect.Lists;
import com.google.common.io.Closeables;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class FullTest {
    private static class TestPosition {
        public final double lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, S12;
        public TestPosition(String line) {
            String[] fields = line.split(" ");
            lat1 = Double.parseDouble(fields[0]);
            lon1 = Double.parseDouble(fields[1]);
            azi1 = Double.parseDouble(fields[2]);
            lat2 = Double.parseDouble(fields[3]);
            lon2 = Double.parseDouble(fields[4]);
            azi2 = Double.parseDouble(fields[5]);
            s12 = Double.parseDouble(fields[6]);
            a12 = Double.parseDouble(fields[7]);
            m12 = Double.parseDouble(fields[8]);
            S12 = Double.parseDouble(fields[9]);
        }
    }
    private static List<TestPosition> data;
    private Geodesic geod;

    @BeforeClass
    public static void readExamples() throws IOException {
        BufferedReader in = null;
        data = Lists.newArrayList();
        try {
            in = new BufferedReader(new InputStreamReader(new GZIPInputStream(FullTest.class.getResourceAsStream("GeodTest-short.dat.gz"))));
            String line;
            while ((line = in.readLine()) != null) {
                data.add(new TestPosition(line));
            }
        } finally {
            Closeables.closeQuietly(in);
        }
    }

    @Before
    public void setGeod() {
        geod = Geodesic.WGS84;
    }

    @Test
    public void testData() {
        assertEquals(10000, data.size());
    }

    @Test
    public void testDirectFromPoint1() {
        Position out = new Position();
        // direct(double lat1, double lon1, double azi1, double s12, int outmask, Position outPosition)
        for (TestPosition p : data) {
            geod.direct(p.lat1, p.lon1, p.azi1, p.s12, Geodesic.Mask.ALL, out);
            assertEquals(p.lat2, out.lat2, 1e-8);
            assertEquals(p.lon2, out.lon2, 1e-8);
            assertEquals(p.azi2, out.azi2, 1e-8);
            assertEquals(p.m12, out.m12, 1e-8);
        }
    }

    @Test
    public void testDirectFromPoint2() {
        Position out = new Position();
        for (TestPosition p : data) {
            geod.direct(p.lat2, p.lon2, p.azi2, -p.s12, Geodesic.Mask.ALL, out);
            assertEquals(p.lat1, out.lat2, 5e-6);
            assertEquals(p.lon1, out.lon2, 5e-6);
            assertEquals(p.azi1, out.azi2, 5e-6);
            assertEquals(-p.m12, out.m12, 5e-6);
        }
    }

    @Test
    public void testInverse() {
        InversePosition out = new InversePosition();
        for (TestPosition p : data) {
            geod.genInverse(p.lat1, p.lon1, p.lat2, p.lon2, Geodesic.Mask.ALL, out);
            // these tolerances seem especially high. are the indicitive of a bug somewhere?
            assertEquals(p.azi1, out.azi1, 1e-3);
            assertEquals(p.azi2, out.azi2, 1e-3);
            assertEquals(p.s12, out.s12, 1e-4);
            assertEquals(p.m12, out.m12, 1e-4);
        }
    }
}
