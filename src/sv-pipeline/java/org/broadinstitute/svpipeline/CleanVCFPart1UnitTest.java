package org.broadinstitute.svpipeline;

public class CleanVCFPart1UnitTest {
    public static void main( final String[] args ) {
        testAsserts();
        testMedianCalculation();
        System.out.println("OK");
    }

    public static void testAsserts() {
        boolean caughtIt = false;
        try {
            assert(false);
        } catch ( final AssertionError ae ) {
            caughtIt = true;
        }
        if ( !caughtIt ) {
            throw new AssertionError("assertions aren't turned on, so you're not testing anything.");
        }
    }

    public static void testMedianCalculation() {
        final int[] counts = new int[4];
        assert(Double.isNaN(CleanVCFPart1.calcMedian(counts)));
        counts[0] = 1;
        assert(CleanVCFPart1.calcMedian(counts) == 0.0);
        counts[1] = 1;
        assert(CleanVCFPart1.calcMedian(counts) == 0.5);
        counts[2] = 1;
        assert(CleanVCFPart1.calcMedian(counts) == 1.0);
        counts[3] = 1;
        assert(CleanVCFPart1.calcMedian(counts) == 1.5);
        counts[2] = 2;
        assert(CleanVCFPart1.calcMedian(counts) == 2.0);
        counts[3] = 4;
        assert(CleanVCFPart1.calcMedian(counts) == 2.5);
        counts[3] = 5;
        assert(CleanVCFPart1.calcMedian(counts) == 3.0);
    }
}
