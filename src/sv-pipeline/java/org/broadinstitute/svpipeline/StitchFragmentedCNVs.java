package org.broadinstitute.svpipeline;

import java.io.*;
import java.util.*;
import org.broadinstitute.svpipeline.VCFParser.*;

/** Read a VCF, and try to stitch together adjacent copy-number variations.
 * Eligible Records (which we call "stitchable") must meet certain criteria as specified by the
 * isStitchable method of the StitchableIterator.
 * If two stitchables overlap appropriately, and all their samples have identical genotypes, we can
 * replace the first one by adding on the interval covered by the second one.
 */
public class StitchFragmentedCNVs {
    // These 3 values will always be overwritten, but are initialized to reasonable defaults as documentation
    private static double PAD_FACTOR = .2;
    static int MAX_PAD = 200000; // visible for testing
    private static double MAX_OVERLAP_FACTOR = .2;

    // A new "end position" for the disposition map signalling that the record is to be removed because it was combined
    // with another record.
    private static final int ENDPOS_REMOVED_RECORD = -1;

    // relevant INFO field keys and values--these are visible for testing
    static final ByteSequence END = new ByteSequence("END");
    static final ByteSequence SVLEN = new ByteSequence("SVLEN");
    static final ByteSequence SVTYPE = new ByteSequence("SVTYPE");
    static final ByteSequence SVTYPE_DEL = new ByteSequence("DEL");
    static final ByteSequence SVTYPE_DUP = new ByteSequence("DUP");
    static final ByteSequence MULTIALLELIC = new ByteSequence("MULTIALLELIC");
    static final ByteSequence EVIDENCE = new ByteSequence("EVIDENCE");
    static final ByteSequence EVIDENCE_RD = new ByteSequence("RD");
    static final ByteSequence EVIDENCE_SR = new ByteSequence("SR");
    static final ByteSequence EVIDENCE_PE = new ByteSequence("PE");
    static final ByteSequence EVIDENCE_BAF = new ByteSequence("BAF");

    public static void main( final String[] args ) {
        if ( args.length != 4 ) {
            System.err.println("Usage: java StitchFragmentedCNVs PAD% MAXPAD OVRLAP% VCFFILE");
            System.err.println("E.g.:  java StitchFragmentedCNVs .2 200000 .2 input.vcf.gz");
            System.err.println("Combines neighboring CNVs with matching genotypes into a larger event.");
            System.err.println("Writes an uncompressed vcf to stdout.");
            System.exit(1);
        }

        initCommandLineArgs(args);

        // a map of IDs onto revised ENDs
        final Map<ByteSequence, Integer> disposition = new HashMap<>(1000);

        // a push-back buffer of previous PaddedIntervals for stitchables that are still in range
        final List<PaddedInterval> intervalList = new LinkedList<>();

        try ( final VCFParser vcfParser = new VCFParser(args[3]) ) {
            while ( vcfParser.hasMetadata() ) {
                vcfParser.nextMetaData();
            }
            ByteSequence currentChromosome = null;
            while ( vcfParser.hasRecord() ) {
                final Record record = vcfParser.nextRecord();
                if ( isStitchable(record) ) {
                    if ( !record.getChromosome().equals(currentChromosome) ) {
                        intervalList.clear();
                        currentChromosome = record.getChromosome();
                    }
                    final PaddedInterval currentInterval = new PaddedInterval(record);
                    final ListIterator<PaddedInterval> previousIntervals = intervalList.listIterator();
                    boolean recordRemoved = false;
                    while ( previousIntervals.hasNext() ) {
                        final PaddedInterval previousInterval = previousIntervals.next();
                        final PaddedInterval revisedInterval;
                        if ( previousInterval.doneStitching(currentInterval) ) {
                            previousIntervals.remove();
                        } else if ( (revisedInterval = previousInterval.stitchTo(currentInterval)) != null ) {
                            previousIntervals.set(revisedInterval);
                            disposition.put(revisedInterval.getRecord().getID(), revisedInterval.getVCFEnd());
                            disposition.put(record.getID(), ENDPOS_REMOVED_RECORD);
                            recordRemoved = true;
                        }
                    }
                    if ( !recordRemoved ) {
                        intervalList.add(currentInterval);
                    }
                }
            }
        }
        intervalList.clear();

        try ( final OutputStream os = new BufferedOutputStream(new FileOutputStream(FileDescriptor.out));
                final VCFParser vcfParser = new VCFParser(args[3])) {
            while ( vcfParser.hasMetadata() ) {
                vcfParser.nextMetaData().write(os);
            }
            while ( vcfParser.hasRecord() ) {
                final Record record = vcfParser.nextRecord();
                final Integer endPosObj = disposition.get(record.getID());
                if ( endPosObj == null ) {
                    record.write(os);
                } else {
                    final int endPos = endPosObj;
                    if ( endPos != ENDPOS_REMOVED_RECORD ) {
                        final InfoField infoField = record.getInfo();
                        infoField.put(END, new ByteSequence(Integer.toString(endPos)));
                        final int svLength = endPos + 1 - record.getPosition();
                        infoField.put(SVLEN, new ByteSequence(Integer.toString(svLength)));
                        record.write(os);
                    }
                }
            }
        } catch ( final IOException ioe ) {
            throw new MalformedVCFException("can't write revised vcf", ioe);
        }
    }

    private static void initCommandLineArgs( final String[] args ) {
        try {
            PAD_FACTOR = Double.parseDouble(args[0]);
        } catch ( final NumberFormatException nfe ) {
            System.err.println("Can't interpret 1st argument (padding fraction) as a floating point number.");
            System.exit(2);
        }
        if ( PAD_FACTOR < 0.0 ) {
            System.err.println("First argument should be a padding fraction >= 0.");
            System.exit(2);
        }
        try {
            MAX_PAD = Integer.parseInt(args[1]);
        } catch ( final NumberFormatException nfe ) {
            System.err.println("Can't interpret 2nd argument (maximum padding in bases) as an integer.");
            System.exit(2);
        }
        if ( MAX_PAD < 0 ) {
            System.err.println("Second argument must be a maximum padding in bases >= 0.");
            System.exit(2);
        }
        try {
            MAX_OVERLAP_FACTOR = Double.parseDouble(args[0]);
        } catch ( final NumberFormatException nfe ) {
            System.err.println("Can't interpret 3rd argument (maximum overlap fraction) as a floating point number.");
            System.exit(2);
        }
        if ( MAX_OVERLAP_FACTOR < 0.0 || MAX_OVERLAP_FACTOR > 1.0 ) {
            System.err.println("Third argument should be a maximum overlap fraction between 0 and 1.");
            System.exit(2);
        }
    }

    // VisibleForTesting
    static boolean isStitchable( final Record record ) {
        final CompoundField filterField = record.getFilter();
        for ( final ByteSequence filter : filterField ) {
            if ( MULTIALLELIC.equals(filter) ) {
                return false;
            }
        }

        final InfoField infoField = record.getInfo();
        final ByteSequence svType = infoField.get(SVTYPE);
        if ( !SVTYPE_DEL.equals(svType) && !SVTYPE_DUP.equals(svType) ) {
            return false;
        }

        // you can't be a stitchable if you don't have an "END" info field.
        // code elsewhere assumes it can grab this value without checking for its existence
        final ByteSequence endValue = infoField.get(END);
        if ( endValue == null || endValue.asInt() == ByteSequence.MISSING_VALUE ) {
            return false;
        }

        final ByteSequence evidence = infoField.get(EVIDENCE);
        if ( evidence == null ) {
            return false;
        }
        return !evidence.contains(EVIDENCE_PE) && !evidence.contains(EVIDENCE_SR) &&
                (evidence.contains(EVIDENCE_RD) || evidence.contains(EVIDENCE_BAF));
    }

    /** A little helper class to do padding and overlap calculations
     *  Note: this class uses half-open intervals, unlike a vcf */
    final static class PaddedInterval { // visible for testing
        private final Record record;
        private final int start;
        private final int end;
        private final int padding;
        private final int maxOverlap;
        private final ByteSequence eventType;

        public PaddedInterval( final Record record ) {
            this.record = record;
            this.start = record.getPosition();
            this.end = record.getInfo().get(END).asInt() + 1;
            final int length = end - start;
            this.padding = Math.min(MAX_PAD, (int)(length * PAD_FACTOR));
            this.maxOverlap = (int)(length * MAX_OVERLAP_FACTOR);
            this.eventType = record.getInfo().get(SVTYPE);
        }

        private PaddedInterval( final PaddedInterval upstream, final PaddedInterval downstream ) {
            this.record = upstream.record;
            this.start = upstream.start;
            this.end = downstream.end;
            this.padding = Math.max(upstream.padding, downstream.padding);
            this.maxOverlap = Math.max(upstream.maxOverlap, downstream.maxOverlap);
            this.eventType = upstream.eventType;
        }

        public int getPaddedStart() { return start - padding; }
        public int getPaddedEnd() { return end + padding; }
        public Record getRecord() { return record; }

        /** Returns true if we're done trying to stitch this interval.  Criterion is that the
         * padded end of this interval more than MAX_PAD bases away from the start of the
         * currentInterval.  So this one is definitely disjoint (regardless of its length), and that
         * will also be true of all subsequent intervals (since they're in sorted order on the
         * starting interval.
         */
        public boolean doneStitching( final PaddedInterval currentInterval ) {
            return getPaddedEnd() < currentInterval.start - MAX_PAD;
        }

        /** Returns an expanded interval if possible, otherwise null. */
        public PaddedInterval stitchTo( final PaddedInterval downstreamInterval ) {
            if ( !eventType.equals(downstreamInterval.eventType) ) {
                return null;
            }

            // Check that the padded intervals overlap.
            // Only have to check one end, because we know the downstream interval starts as late
            // or later than this one.
            if ( getPaddedEnd() <= downstreamInterval.getPaddedStart() ) {
                return null;
            }

            // But the unpadded intervals mustn't overlap too much.
            // Note that the calculated overlap can be negative (they don't actually overlap),
            //   but that's OK.
            final int overlap = Math.min(end, downstreamInterval.end) - downstreamInterval.start;
            if ( overlap > maxOverlap || overlap > downstreamInterval.maxOverlap ) {
                return null;
            }

            if ( !genotypesMatch(record.getGenotypes(), downstreamInterval.record.getGenotypes()) ) {
                return null;
            }

            return new PaddedInterval(this, downstreamInterval);
        }

        public int getVCFEnd() { return end - 1; }

        private static boolean genotypesMatch( final List<CompoundField> genotypes1,
                                               final List<CompoundField> genotypes2 ) {
            final int nGTs = genotypes1.size();
            if ( genotypes2.size() != nGTs ) {
                throw new IllegalStateException("records have a different number of genotypes");
            }
            for ( int idx = 0; idx != nGTs; ++idx ) {
                final ByteIterator itr1 = genotypes1.get(idx).getValue().iterator();
                final ByteIterator itr2 = genotypes2.get(idx).getValue().iterator();
                byte b1;
                do {
                    b1 = itr1.hasNext() ? itr1.next() : (byte)':';
                    final byte b2 = itr2.hasNext() ? itr2.next() : (byte)':';
                    if ( b1 != b2 ) return false;
                } while ( b1 != ':' );
            }
            return true;
        }
    }
}
