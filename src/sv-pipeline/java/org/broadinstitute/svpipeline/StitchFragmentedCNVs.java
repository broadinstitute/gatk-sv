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
    private static int MAX_PAD = 200000;
    private static double MAX_OVERLAP_FACTOR = .2;

    // A new "end position" for the disposition map signalling that the record is to be removed because it was combined
    // with another record.
    private static final int ENDPOS_REMOVED_RECORD = -1;

    private static final ByteSequence END = new ByteSequence("END");
    private static final ByteSequence SVLEN = new ByteSequence("SVLEN");
    private static final ByteSequence SVTYPE = new ByteSequence("SVTYPE");

    public static void main( final String[] args ) {
        if ( args.length != 4 ) {
            System.err.println("Usage: java StitchFragmentedCNVs PAD% MAXPAD OVRLAP% VCFFILE");
            System.err.println("E.g.:  java StitchFragmentedCNVs .2 200000 .2 input.vcf.gz");
            System.err.println("Combines neighboring CNVs with matching genotypes into a larger event.");
            System.err.println("Writes an uncompressed vcf to stdout.");
            System.exit(1);
        }

        initCommandLineArgs(args);

        final Map<ByteSequence, Integer> disposition = new HashMap<>(1000);
        try ( final VCFParser vcfParser = new VCFParser(args[3]) ) {
            while ( vcfParser.hasMetadata() ) {
                vcfParser.nextMetaData();
            }
            final StitchableIterator sItr = new StitchableIterator(vcfParser);
            Record stitchableRecord;
            while ( (stitchableRecord = sItr.nextSubject()) != null ) {
                findExtension(stitchableRecord, sItr, disposition);
            }
        } catch ( final IOException ioe ) {
            throw new MalformedVCFException("can't read vcf", ioe);
        }

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

    /** Look for a stitchable downstream of the subject that can be joined to it to
     *  make a larger event. */
    private static void findExtension( final Record stitchable,
                                       final StitchableIterator sItr,
                                       final Map<ByteSequence, Integer> disposition )
            throws IOException {

        final PaddedInterval originalPaddedInterval = new PaddedInterval(stitchable);
        PaddedInterval paddedInterval = originalPaddedInterval;

        // sItr.hasNext returns false at EOF, or when the next record is too far away to
        // overlap the subject
        while ( sItr.hasNext() ) {
            final Record record = sItr.next();
            final PaddedInterval paddedInterval2 = new PaddedInterval(record);
            if ( paddedInterval.canCoalesceWith(paddedInterval2) &&
                    genotypesMatch(stitchable.getGenotypes(), record.getGenotypes()) ) {
                paddedInterval = paddedInterval2;
                disposition.put(record.getID(), ENDPOS_REMOVED_RECORD);
            }
        }

        if ( paddedInterval != originalPaddedInterval ) {
            final int endPos = paddedInterval.getVCFEnd();
            disposition.put(stitchable.getID(), endPos);
        }
    }

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

    /** A little helper class to do padding and overlap calculations
     *  Note: this class uses half-open intervals, unlike a vcf */
    private final static class PaddedInterval {
        private final int start;
        private final int end;
        private final int padding;
        private final int maxOverlap;
        private final ByteSequence eventType;

        public PaddedInterval( final Record record ) {
            this.start = record.getPosition();
            // getInfoField can't return null -- it's been checked in isStitchable
            // + 1 because vcf has closed intervals, we use half-open
            this.end = record.getInfo().get(END).asInt() + 1;
            final int length = end - start;
            this.padding = Math.min(MAX_PAD, (int)(length * PAD_FACTOR));
            this.maxOverlap = (int)(length * MAX_OVERLAP_FACTOR);
            this.eventType = record.getInfo().get(SVTYPE);
        }

        public boolean canCoalesceWith( final PaddedInterval downstreamInterval ) {
            if ( !eventType.equals(downstreamInterval.eventType) ) {
                return false;
            }

            // Check that the padded intervals overlap.
            // Only have to check one end, because we know the downstream interval starts as late
            // or later than this one.
            if ( end + padding <= downstreamInterval.start - downstreamInterval.padding ) {
                return false;
            }

            // but the unpadded intervals mustn't overlap too much
            final int overlap = Math.min(end, downstreamInterval.end) - downstreamInterval.start;
            return overlap < maxOverlap && overlap < downstreamInterval.maxOverlap;
        }

        public int getVCFEnd() { return end - 1; }
    }

    /** Maintains a list of stitchables.
     * This is a kind of double iterator:  At the outer level, calling nextSubject repeatedly
     * until it returns null lets you simply iterate over each stitchable record in the input file.
     * Each time you do so, the inner iterator (hasNext/next) gets reset to iterate over the
     * stitchable records downstream of the subject.  The inner iterator is smart enough to quit
     * (i.e., hasNext will return false) when we've read so far ahead that we can't possibly find
     * a stitchable that can be joined to the subject.
     */
    private final static class StitchableIterator implements Iterator<Record> {
        private final VCFParser vcfParser;
        private final LinkedList<Record> records; // read-ahead for stitchables downstream of subject
        private Record nextRecord; // next record for inner iteration
        private Iterator<Record> innerIterator;
        private ByteSequence subjectChromosome;
        private int subjectMinNoOverlapPosition; // far enough downstream that MAX_PAD will ensure there's no overlap

        private static final ByteSequence MULTIALLELIC = new ByteSequence("MULTIALLELIC");
        private static final ByteSequence SVTYPE_DEL = new ByteSequence("DEL");
        private static final ByteSequence SVTYPE_DUP = new ByteSequence("DUP");
        private static final ByteSequence EVIDENCE = new ByteSequence("EVIDENCE");
        private static final ByteSequence EVIDENCE_RD = new ByteSequence("RD");
        private static final ByteSequence EVIDENCE_SR = new ByteSequence("SR");
        private static final ByteSequence EVIDENCE_PE = new ByteSequence("PE");
        private static final ByteSequence EVIDENCE_BAF = new ByteSequence("BAF");

        public StitchableIterator( final VCFParser vcfParser ) {
            this.vcfParser = vcfParser;
            this.records = new LinkedList<>();
        }

        /** return the next stitchable subject */
        public Record nextSubject() throws IOException {
            if ( !records.isEmpty() ) {
                final Record next = records.removeFirst();
                innerIterator = records.iterator();
                return setSubject(next);
            }

            while ( vcfParser.hasRecord() ) {
                final Record next = vcfParser.nextRecord();
                if ( isStitchable(next) ) {
                    innerIterator = Collections.emptyIterator();
                    return setSubject(next);
                }
            }

            return null;
        }

        /** Is there another stitchable downstream of the subject that is within joining range? */
        @Override public boolean hasNext() {
            if ( nextRecord != null ) {
                return true;
            }

            if ( innerIterator == null ) {
                return false;
            }

            if ( innerIterator.hasNext() ) {
                return setupInRange(innerIterator.next());
            }

            innerIterator = Collections.emptyIterator();
            while ( vcfParser.hasRecord() ) {
                final Record record = vcfParser.nextRecord();
                if ( isStitchable(record) ) {
                    records.addLast(record);
                    return setupInRange(record);
                }
            }
            return false;
        }

        @Override public Record next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException();
            }
            final Record result = nextRecord;
            nextRecord = null;
            return result;
        }

        private boolean setupInRange( final Record record ) {
            if ( record.getChromosome().equals(subjectChromosome) &&
                    record.getPosition() < subjectMinNoOverlapPosition ) {
                nextRecord = record;
                return true;
            }
            nextRecord = null;
            innerIterator = null;
            return false;
        }

        private static boolean isStitchable( final Record record ) {
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

        private Record setSubject( final Record record ) {
            subjectChromosome = record.getChromosome();
            final int start = record.getPosition();
            final int end = record.getInfo().get(END).asInt();
            subjectMinNoOverlapPosition =
                    end + Math.min(MAX_PAD, (int)(PAD_FACTOR * (end - start))) + MAX_PAD;
            return record;
        }
    }
}
