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
    private static final ByteSequence END = new ByteSequence("END");
    private static final ByteSequence SVLEN = new ByteSequence("SVLEN");

    // These 3 values will always be overwritten, but are initialized to reasonable defaults as documentation
    private static double PAD_FACTOR = .2;
    private static int MAX_PAD = 200000;
    private static double MAX_OVERLAP_FACTOR = .2;

    public static void main( final String[] args ) {
        if ( args.length != 4 ) {
            System.err.println("Usage: java StitchFragmentedCNVs PAD% MAXPAD OVRLAP% VCFFILE");
            System.err.println("E.g.:  java StitchFragmentedCNVs .2 200000 .2 input.vcf.gz");
            System.err.println("Combines neighboring CNVs with matching genotypes into a larger event.");
            System.err.println("Writes an uncompressed vcf to stdout.");
            System.exit(1);
        }

        initCommandLineArgs(args);

        try ( final OutputStream os
                      = new BufferedOutputStream(new FileOutputStream(FileDescriptor.out)) ) {
            try ( final VCFParser vcfParser = new VCFParser(args[3]) ) {
                while ( vcfParser.hasMetadata() ) {
                    vcfParser.nextMetaData().write(os);
                }
                final StitchableIterator sItr = new StitchableIterator(vcfParser);
                Record stitchableRecord;
                while ( (stitchableRecord = sItr.nextSubject(os)) != null ) {
                    findExtension(stitchableRecord, sItr);
                    stitchableRecord.write(os);
                }
            }
        } catch ( final IOException ioe ) {
            throw new MalformedVCFException("can't write to stdout", ioe);
        }
    }

    /** Look for a stitchable downstream of the subject that can be joined to it to
     *  make a larger event. */
    private static void findExtension( final Record stitchable,
                                       final StitchableIterator sItr ) throws IOException {
        final PaddedInterval originalPaddedInterval = new PaddedInterval(stitchable);
        PaddedInterval paddedInterval = originalPaddedInterval;

        // sItr.hasNext returns false at EOF, or when the next record is too far away to
        // overlap the subject
        while ( sItr.hasNext() ) {
            final Record record = sItr.next();
            final PaddedInterval paddedInterval2 = new PaddedInterval(record);
            if ( paddedInterval.canCoalesceWith(paddedInterval2) &&
                    stitchable.getGenotypes().equals(record.getGenotypes()) ) {
                paddedInterval = paddedInterval2;
                sItr.remove();
            }
        }

        if ( paddedInterval != originalPaddedInterval ) {
            final int endPos = paddedInterval.getVCFEnd();
            final InfoField infos = stitchable.getInfo();
            infos.put(END, new ByteSequence(Integer.toString(endPos)));
            final int svLength = endPos + 1 - stitchable.getPosition();
            infos.put(SVLEN, new ByteSequence(Integer.toString(svLength)));
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

    /** A little helper class to do padding and overlap calculations
     *  Note: this class uses half-open intervals, unlike a vcf */
    private final static class PaddedInterval {
        private final int start;
        private final int end;
        private final int padding;
        private final int maxOverlap;

        public PaddedInterval( final Record record ) {
            this.start = record.getPosition();
            // getInfoField can't return null -- it's been checked in isStitchable
            // + 1 because vcf has closed intervals, we use half-open
            this.end = record.getInfo().get(END).asInt() + 1;
            final int length = end - start;
            this.padding = Math.min(MAX_PAD, (int)(length * PAD_FACTOR));
            this.maxOverlap = (int)(length * MAX_OVERLAP_FACTOR);
        }

        public boolean canCoalesceWith( final PaddedInterval downstreamInterval ) {
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

    /** As we go through the VCF we create a new Chunk whenever we encounter a stitchable record.
     *  So, a Chunk consists of a mess of non-stitchables, and a trailing stitchable. */
    private final static class Chunk {
        private final List<Record> nonStitchables;
        private final Record stitchable;

        public Chunk( final List<Record> nonStitchables,
                      final Record stitchable ) {
            this.nonStitchables = nonStitchables;
            this.stitchable = stitchable;
        }

        public List<Record> getNonStitchables() { return nonStitchables; }
        public Record getStitchable() { return stitchable; }
    }

    /** Maintains a list of chunks so that the client just sees the stitchables, while making sure
     * that the record ordering of the input file is maintained.
     * This is a kind of double iterator:  At the outer level, calling nextSubject repeatedly
     * until it returns null lets you simply iterate over each stitchable record in the input file.
     * Each time you do so, the inner iterator (hasNext/next) gets reset to iterate over the
     * stitchable records downstream of the subject.  The inner iterator is smart enough to quit
     * (i.e., hasNext will return false) when we've read so far ahead that we can't possibly find
     * a stitchable that can be joined to the subject.
     */
    private final static class StitchableIterator implements Iterator<Record> {
        private final VCFParser vcfParser;
        private final List<Chunk> chunks;
        private int subjectIndex;
        private ByteSequence subjectChromosome;
        private int subjectMinNoOverlapPosition; // far enough downstream that MAX_PAD will ensure there's no overlap
        private int iterationIndex;
        private Record nextRecord; // this is a pushback for a record that's too far downstream

        private static final ByteSequence MULTIALLELIC = new ByteSequence("MULTIALLELIC");
        private static final ByteSequence SVTYPE = new ByteSequence("SVTYPE");
        private static final ByteSequence SVTYPE_DEL = new ByteSequence("DEL");
        private static final ByteSequence SVTYPE_DUP = new ByteSequence("DUP");
        private static final ByteSequence EVIDENCE = new ByteSequence("EVIDENCE");
        private static final ByteSequence EVIDENCE_RD = new ByteSequence("RD");
        private static final ByteSequence EVIDENCE_SR = new ByteSequence("SR");
        private static final ByteSequence EVIDENCE_PE = new ByteSequence("PE");
        private static final ByteSequence EVIDENCE_BAF = new ByteSequence("BAF");

        public StitchableIterator( final VCFParser vcfParser ) {
            this.vcfParser = vcfParser;
            this.chunks = new ArrayList<>();
        }

        /** write the non-stitchables that precede the first stitchable,
         *  and return the next stitchable */
        public Record nextSubject( final OutputStream os ) throws IOException {
            final int nChunks = chunks.size();
            while ( subjectIndex < nChunks ) {
                final Chunk chunk = chunks.get(subjectIndex);
                // clean out the chunks as we use them:  in coordinate-dense vcfs the chunks array
                // can get quite large, and, especially in vcfs with lots of sample, a chunk can
                // occupy quite a large amount of memory.  we want to release the chunks for garbage
                // collection ASAP to control memory use.
                chunks.set(subjectIndex, null);
                iterationIndex = ++subjectIndex;
                for ( final Record rec : chunk.getNonStitchables() ) {
                    rec.write(os);
                }
                final Record stitchable = chunk.getStitchable();
                if ( stitchable != null ) {
                    return setSubject(stitchable);
                }
            }

            // there are no more chunks to serve as subjects, reset the queue
            chunks.clear();
            subjectIndex = iterationIndex = 0;

            while ( nextRecord != null || vcfParser.hasRecord() ) {
                final Record record = nextRecord != null ? nextRecord: vcfParser.nextRecord();
                nextRecord = null;
                if ( isStitchable(record) ) {
                    return setSubject(record);
                }
                record.write(os);
            }
            return null;
        }

        /** Is there another stitchable downstream of the subject that is within joining range? **/
        @Override public boolean hasNext() {
            final int nChunks = chunks.size();
            while ( iterationIndex < nChunks ) {
                final Record stitchable = chunks.get(iterationIndex).getStitchable();
                if ( stitchable != null ) {
                    return true;
                }
                ++iterationIndex;
            }
            if ( nextRecord != null || vcfParser.hasRecord() ) {
                List<Record> nonStitchables = null;
                do {
                    final Record record =
                            nextRecord != null ? nextRecord : vcfParser.nextRecord();
                    nextRecord = null;
                    if ( !record.getChromosome().equals(subjectChromosome) ||
                            record.getPosition() >= subjectMinNoOverlapPosition ) {
                        nextRecord = record;
                        if ( nonStitchables != null ) {
                            chunks.add(new Chunk(nonStitchables, null));
                        }
                        return false;
                    }
                    if ( isStitchable(record) ) {
                        if ( nonStitchables == null ) {
                            nonStitchables = Collections.emptyList();
                        }
                        chunks.add(new Chunk(nonStitchables, record));
                        return true;
                    }
                    if ( nonStitchables == null ) {
                        nonStitchables = new ArrayList<>();
                    }
                    nonStitchables.add(record);
                } while ( vcfParser.hasRecord() );
                chunks.add(new Chunk(nonStitchables, null));
            }
            return false;
        }

        @Override public Record next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException();
            }
            return chunks.get(iterationIndex++).getStitchable();
        }

        @Override public void remove() {
            final int idx = iterationIndex - 1;
            chunks.set(idx, new Chunk(chunks.get(idx).getNonStitchables(), null));
        }

        private static boolean isStitchable( final Record record ) {
            if ( MULTIALLELIC.equals(record.getFilter().getValue()) ) return false;
            final InfoField info = record.getInfo();
            final ByteSequence svType = info.get(SVTYPE);
            if ( !SVTYPE_DEL.equals(svType) && !SVTYPE_DUP.equals(svType) ) return false;
            // you can't be a stitchable if you don't have an "END" info field.
            // code elsewhere assumes it can grab this value without checking for its existence
            final ByteSequence end = info.get(END);
            if ( end == null || end.asInt() == ByteSequence.MISSING_VALUE ) return false;
            final ByteSequence evidence = info.get(EVIDENCE);
            if ( evidence == null ) return false;
            return !evidence.contains(EVIDENCE_PE) && !evidence.contains(EVIDENCE_SR) &&
                    (evidence.contains(EVIDENCE_RD) || evidence.contains(EVIDENCE_BAF));
        }

        private Record setSubject( final Record record ) {
            subjectChromosome = record.getChromosome();
            final int start = record.getPosition();
            final int end = record.getInfo().get(END).asInt(); // can't be null -- checked in isStitchable
            subjectMinNoOverlapPosition =
                    end + Math.min(MAX_PAD, (int)(PAD_FACTOR * (end - start))) + MAX_PAD;
            return record;
        }
    }
}
