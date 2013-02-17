package com.github.mlc.jgeographic;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;
import com.google.common.base.CharMatcher;
import com.google.common.base.Function;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.io.Closeables;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class FullTest {
    private static class TestPosition {
        private static final Splitter SPACE_SPLITTER = Splitter.on(CharMatcher.BREAKING_WHITESPACE).omitEmptyStrings();
        public final double lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, S12;
        public TestPosition(String line) {
            Iterator<String> strs = SPACE_SPLITTER.split(line).iterator();
            Iterator<Double> doubles = Iterators.transform(strs, new Function<String, Double>() {
                @Override
                public Double apply(String input) {
                    return Double.parseDouble(input);
                }
            });
            lat1 = doubles.next();
            lon1 = doubles.next();
            azi1 = doubles.next();
            lat2 = doubles.next();
            lon2 = doubles.next();
            azi2 = doubles.next();
            s12 = doubles.next();
            a12 = doubles.next();
            m12 = doubles.next();
            S12 = doubles.next();
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
