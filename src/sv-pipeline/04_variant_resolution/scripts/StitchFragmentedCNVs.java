import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/** Read a VCF, and try to stitch together adjacent copy-number variations.
 * Eligible Records (which we call "stitchable") must meet certain criteria as specified by the
 * isStitchable method of the StitchableIterator.
 * If two stitchables overlap in a weird way that seems wrong to me (see line 54), and all their
 * samples have identical genotypes, we can replace the first one by adding on the interval covered
 * by the second one.
 */
public class StitchFragmentedCNVs {
    private static final VCFParser.ByteSequence END = new VCFParser.ByteSequence("END");
    private static final VCFParser.ByteSequence SVLEN = new VCFParser.ByteSequence("SVLEN");
    private static final double PAD_FACTOR = .2;
    private static final int MAX_PAD = 200000;
    private static final String GZ = ".gz";

    public static void main( final String[] args ) {
        try ( final OutputStream os = createOutputStream(args[1]) ) {
            try ( final VCFParser vcfParser = new VCFParser(args[0]) ) {
                while ( vcfParser.hasMetadata() ) {
                    vcfParser.nextMetaData().write(os);
                }
                final StitchableIterator sItr = new StitchableIterator(vcfParser);
                VCFParser.Record stitchableRecord;
                while ( (stitchableRecord = sItr.nextSubject(os)) != null ) {
                    findExtension(stitchableRecord, sItr);
                    stitchableRecord.write(os);
                }
            }
        } catch ( final IOException ioe ) {
            throw new VCFParser.MalformedVCFException("can't write output to " + args[1], ioe);
        }
    }

    private static void findExtension( final VCFParser.Record stitchable,
                                       final StitchableIterator sItr ) throws IOException {
        final int startPos = stitchable.getPosition();
        final VCFParser.ByteSequence endField = stitchable.getInfoField(END);
        final int originalEndPos = endField.asInt(); // can't be null -- checked in isStitchable
        int endPos = originalEndPos;
        int pad = Math.min(MAX_PAD, (int)(PAD_FACTOR * (endPos - startPos)));
        final int padStartPos = Math.max(1, startPos - pad);
        int padEndPos = endPos + pad;

        // sItr.hasNext returns false at EOF, or when the next record is too far away to overlap the subject
        while ( sItr.hasNext() ) {
            final VCFParser.Record record = sItr.next();
            final int startPos2 = record.getPosition();
            final int endPos2 = record.getInfoField(END).asInt(); // can't be null -- checked in isStitchable
            final int pad2 = Math.min(MAX_PAD, (int)(PAD_FACTOR * (endPos2 - startPos2)));
            final int padStartPos2 = Math.max(1, startPos2 - pad2);
            final int padEndPos2 = endPos2 + pad2;
            final int overlap = endPos - startPos2; // this seems clearly wrong. need clarification.
            if ( overlap < pad && overlap < pad2 &&
                    Math.max(padStartPos, padStartPos2) < Math.min(padEndPos, padEndPos2) &&
                    genotypesMatch(stitchable, record) ) {
                endPos = Math.max(endPos, endPos2);
                pad = Math.min(MAX_PAD, (int)(PAD_FACTOR * (endPos - startPos)));
                padEndPos = endPos + pad;
                sItr.remove();
            }
        }

        if ( endPos != originalEndPos ) {
            stitchable.setInfoField(endField, new VCFParser.ByteSequence(Integer.toString(endPos)));
            final VCFParser.ByteSequence svLenField = stitchable.getInfoField(SVLEN);
            if ( svLenField != null ) {
                final VCFParser.ByteSequence svLenValue =
                        new VCFParser.ByteSequence(Integer.toString(endPos - startPos));
                stitchable.setInfoField(svLenField, svLenValue);
            }
        }
    }

    private static boolean genotypesMatch( final VCFParser.Record rec1, final VCFParser.Record rec2 ) {
        final List<VCFParser.ByteSequence> gt1 = rec1.getGenotypes();
        final List<VCFParser.ByteSequence> gt2 = rec2.getGenotypes();
        final int nGTs = gt1.size();
        if ( gt2.size() != nGTs ) {
            throw new IllegalStateException("two records have a different number of genotypes");
        }
        for ( int idx = 0; idx != nGTs; ++idx ) {
            final VCFParser.ByteIterator itr1 = gt1.get(idx).iterator();
            final VCFParser.ByteIterator itr2 = gt2.get(idx).iterator();
            byte b1;
            do {
                b1 = itr1.hasNext() ? itr1.next() : -1;
                final byte b2 = itr2.hasNext() ? itr2.next() : -1;
                if ( b1 != b2 ) return false;
            } while ( b1 != ':' );
        }
        return true;
    }

    private static OutputStream createOutputStream( final String pathName ) throws IOException {
        final OutputStream os = pathName.equals("-") ? System.out : new FileOutputStream(pathName);
        return new BufferedOutputStream(pathName.endsWith(GZ) ? new GZIPOutputStream(os) : os);
    }

    /** As we go through the VCF we create a new Chunk whenever we encounter a stitchable record.
     * So, a Chunk consists of a mess of non-stitchables, and a trailing stitchable. */
    private final static class Chunk {
        private final List<VCFParser.Record> nonStitchables;
        private final VCFParser.Record stitchable;

        public Chunk( final List<VCFParser.Record> nonStitchables, final VCFParser.Record stitchable ) {
            this.nonStitchables = nonStitchables;
            this.stitchable = stitchable;
        }

        public List<VCFParser.Record> getNonStitchables() { return nonStitchables; }
        public VCFParser.Record getStitchable() { return stitchable; }
    }

    /** maintains a list of chunks so that the client just sees the stitchables.
     * also acts like a queue, in that you can remove the first Chunk by calling nextSubject.
     */
    private final static class StitchableIterator implements Iterator<VCFParser.Record> {
        private final VCFParser vcfParser;
        private final List<Chunk> chunks;
        private int subjectIndex;
        private VCFParser.ByteSequence subjectChromosome;
        private int subjectMinNoOverlapPosition; // far enough downstream that MAX_PAD will ensure there's no overlap
        private int iterationIndex;

        private static final VCFParser.ByteSequence MULTIALLELIC = new VCFParser.ByteSequence("MULTIALLELIC");
        private static final VCFParser.ByteSequence SVTYPE = new VCFParser.ByteSequence("SVTYPE");
        private static final VCFParser.ByteSequence SVTYPE_DEL = new VCFParser.ByteSequence("DEL");
        private static final VCFParser.ByteSequence SVTYPE_DUP = new VCFParser.ByteSequence("DUP");
        private static final VCFParser.ByteSequence EVIDENCE = new VCFParser.ByteSequence("EVIDENCE");
        private static final String EVIDENCE_RD = "RD";
        private static final String EVIDENCE_SR = "SR";
        private static final String EVIDENCE_PE = "PE";
        private static final String EVIDENCE_BAF = "BAF";

        public StitchableIterator( final VCFParser vcfParser ) {
            this.vcfParser = vcfParser;
            this.chunks = new ArrayList<>();
            this.subjectIndex = 0;
            this.iterationIndex = 0;
        }

        /** write the non-stitchables that preceed the first stitchable, and return the stitchable */
        public VCFParser.Record nextSubject( final OutputStream os ) throws IOException {
            final int nChunks = chunks.size();
            while ( subjectIndex < nChunks ) {
                final Chunk chunk = chunks.get(subjectIndex);
                chunks.set(subjectIndex++, null);
                iterationIndex = subjectIndex + 1;
                for ( final VCFParser.Record rec : chunk.getNonStitchables() ) {
                    rec.write(os);
                }
                final VCFParser.Record stitchable = chunk.getStitchable();
                if ( stitchable != null ) {
                    return setSubject(stitchable);
                }
            }
            // there are no more chunks to serve as subjects, reset the queue
            chunks.clear();
            subjectIndex = iterationIndex = 0;

            while ( vcfParser.hasRecord() ) {
                final VCFParser.Record record = vcfParser.nextRecord();
                if ( isStitchable(record) ) {
                    return setSubject(record);
                }
                record.write(os);
            }
            return null;
        }

        @Override public boolean hasNext() {
            final int nChunks = chunks.size();
            while ( iterationIndex < nChunks ) {
                final VCFParser.Record stitchable = chunks.get(iterationIndex).getStitchable();
                if ( stitchable != null ) {
                    return true;
                }
                ++iterationIndex;
            }
            if ( vcfParser.hasRecord() ) {
                final List<VCFParser.Record> nonStitchables = new ArrayList<>();
                do {
                    final VCFParser.Record record = vcfParser.nextRecord();
                    if ( isStitchable(record) ) {
                        chunks.add(new Chunk(nonStitchables, record));
                        return true;
                    }
                    nonStitchables.add(record);
                    // don't just build up Chunks indefinitely:
                    // break when we've read ahead too far downstream of the subject
                    if ( !record.getChromosome().equals(subjectChromosome) ||
                            record.getPosition() >= subjectMinNoOverlapPosition ) {
                        break;
                    }
                } while ( vcfParser.hasRecord() );
                chunks.add(new Chunk(nonStitchables, null));
            }
            return false;
        }

        @Override public VCFParser.Record next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException();
            }
            return chunks.get(iterationIndex++).getStitchable();
        }

        @Override public void remove() {
            final int idx = iterationIndex - 1;
            chunks.set(idx, new Chunk(chunks.get(idx).getNonStitchables(), null));
        }

        private static boolean isStitchable( final VCFParser.Record record ) {
            if ( MULTIALLELIC.equals(record.getFilter()) ) return false;
            final Map<VCFParser.ByteSequence, VCFParser.ByteSequence> infoMap = record.getInfoAsMap();
            final VCFParser.ByteSequence svType = infoMap.get(SVTYPE);
            if ( !SVTYPE_DEL.equals(svType) && !SVTYPE_DUP.equals(svType) ) return false;
            final VCFParser.ByteSequence end = infoMap.get(END);
            if ( end == null || end.asInt() == VCFParser.ByteSequence.MISSING_VALUE ) return false;
            final VCFParser.ByteSequence evidence = infoMap.get(EVIDENCE);
            if ( evidence == null ) return false;
            final String evStr = evidence.toString();
            return !evStr.contains(EVIDENCE_PE) && !evStr.contains(EVIDENCE_SR) &&
                    (evStr.contains(EVIDENCE_RD) || evStr.contains(EVIDENCE_BAF));
        }

        private VCFParser.Record setSubject( final VCFParser.Record record ) {
            subjectChromosome = record.getChromosome();
            final int start = record.getPosition();
            final int end = record.getInfoField(END).asInt(); // can't be null -- checked in isStitchable
            subjectMinNoOverlapPosition =
                    end + Math.min(MAX_PAD, (int)(PAD_FACTOR * (end - start))) + MAX_PAD;
            return record;
        }
    }
}
