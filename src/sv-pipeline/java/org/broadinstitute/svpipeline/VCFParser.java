package org.broadinstitute.svpipeline;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/** Class to parse VCF files.
 * Super simple -- no validation.  Just tokenizes.
 * Aims to avoid parsing fields (like INFO) until necessary.
 * Keeps fields in the form of ByteSequences (like strings, except bytes, not characters).
 *
 * To use:  first iterate over the metadata using hasMetadata() and nextMetaData(), then
 * iterate over the data using hasRecord() and nextRecord().
 */
public class VCFParser implements Closeable {
    private final String pathName;
    private final InputStream is;
    private ByteIterator bufferIterator;

    private static final int BUFFER_SIZE = 64*1024;
    private static final ByteSequence EMPTY_SEQUENCE = new ByteSequence(new byte[0], 0, 0);
    private static final String GZ = ".gz";

    public VCFParser( final String pathName ) {
        if ( pathName == null || "-".equals(pathName) ) {
            this.pathName = "stdin";
            this.is = new BufferedInputStream(new FileInputStream(FileDescriptor.in), BUFFER_SIZE);
        } else {
            this.pathName = pathName;
            try {
                final BufferedInputStream bis =
                        new BufferedInputStream(new FileInputStream(pathName), BUFFER_SIZE);
                this.is = pathName.endsWith(GZ) ? new GZIPInputStream(bis, BUFFER_SIZE) : bis;
            } catch ( final IOException ioe ) {
                throw new MalformedVCFException("can't open " + pathName, ioe);
            }
        }
        if ( !readBuffer() ) {
            throw new MalformedVCFException(this.pathName + " is empty");
        }
    }

    public VCFParser( final InputStream is ) {
        this.pathName = "input VCF";
        this.is = is;
        if ( !readBuffer() ) {
            throw new MalformedVCFException("input VCF is empty");
        }
    }

    public void close() {
        try {
            is.close();
        } catch ( final IOException ioe ) {
            throw new MalformedVCFException("can't close " + pathName, ioe);
        }
    }

    /** process lines beginning with '#' as metadata */
    public boolean hasMetadata() {
        if ( !bufferIterator.hasNext() ) {
            if ( !readBuffer() ) {
                return false;
            }
        }
        return bufferIterator.peek() == '#';
    }

    public Metadata nextMetaData() {
        if ( !hasMetadata() ) {
            throw new NoSuchElementException();
        }
        bufferIterator.skip();
        needData();
        // are we about to parse the column header?
        // it's the only metadata line that doesn't start with "##" but goes:
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT   (sample names)
        if ( bufferIterator.peek() != '#' ) {
            return new ColumnHeaderMetadata(captureColumns());
        }
        bufferIterator.skip();
        // get the key part of the metadata, e.g., INFO in ##INFO=, or contig in ##contig=
        final ByteSequence key = capture('=');
        needData();
        // metadata is just ##key=value, but some keys (like INFO) have multiple attributes
        // as a value.  check to see if we need to tokenize multiple values.
        // multiple values are surrounded by '<' and '>'
        if ( bufferIterator.peek() != '<' ) {
            // nope, simple value.  just grab the rest of the line
            final ByteSequence value = capture('\n');
            return new KeyValueMetadata(key, value);
        }
        bufferIterator.skip();
        // yup.  multiple values. tokenize them.
        return new KeyAttributesMetadata(key, captureAttributes());
    }

    /** once we've had hasMetadata return false (line doesn't start with '#')
     * we're ready to read records.
     * if there are more lines to read, they're all data records. */
    public boolean hasRecord() {
        if ( !bufferIterator.hasNext() ) {
            return readBuffer();
        }
        return true;
    }

    /** tokenize the next line as columns, and stuff them into a Record */
    public Record nextRecord() {
        if ( !hasRecord() ) {
            throw new NoSuchElementException();
        }
        if ( bufferIterator.peek() == '#' ) {
            throw new IllegalStateException("can't read records until metadata is exhausted");
        }
        return new Record(captureColumns());
    }

    /** make sure the bufferIterator has more data, reading another buffer-ful if necessary */
    private void needData() {
        if ( !bufferIterator.hasNext() ) {
            if ( !readBuffer() ) {
                throw new MalformedVCFException("final line truncated in " + pathName);
            }
        }
    }


    private boolean readBuffer() {
        final byte[] buffer = new byte[BUFFER_SIZE];
        try {
            bufferIterator = new ByteIterator(buffer, 0, Math.max(0, is.read(buffer)));
        } catch ( final IOException ioe ) {
            throw new MalformedVCFException("can't read " + pathName);
        }
        return bufferIterator.hasNext();
    }

    /** grab the sequence of bytes up to the specified delimiter */
    private ByteSequence capture( final char delim ) {
        ByteSequence prefix = null;
        int start = bufferIterator.mark();
        do {
            if ( !bufferIterator.hasNext() ) {
                final ByteSequence bs = bufferIterator.getSequence(start);
                if ( bs != null ) {
                    prefix = prefix == null ? bs : new ByteSequence(prefix, bs);
                }
                needData();
                start = 0;
            }
        } while ( bufferIterator.next() != delim );
        final ByteSequence bs = bufferIterator.getSequenceNoDelim(start);
        return prefix == null ?
                (bs == null ? EMPTY_SEQUENCE : bs) :
                (bs == null ? prefix : new ByteSequence(prefix, bs));
    }

    /** tokenize the key/value pairs in a metadata line with multiple attributes
     * ##KEY=<key1=val1,key2=val2,...> */
    private List<KeyValue> captureAttributes() {
        final List<KeyValue> attributes = new ArrayList<>();
        byte finalByte;
        do {
            final ByteSequence key = capture('=');
            needData();
            ByteSequence prefix = null;
            int start = bufferIterator.mark();
            boolean inQuote = false;
            do {
                if ( !bufferIterator.hasNext() ) {
                    final ByteSequence bs = bufferIterator.getSequence(start);
                    if ( bs != null ) {
                        prefix = prefix == null ? bs : new ByteSequence(prefix, bs);
                    }
                    needData();
                    start = 0;
                }
                finalByte = bufferIterator.next();
                if ( finalByte == '"' ) {
                    inQuote = !inQuote;
                }
            } while ( inQuote || (finalByte != ',' && finalByte != '>') );
            final ByteSequence bs = bufferIterator.getSequenceNoDelim(start);
            final ByteSequence value = prefix == null ?
                    (bs == null ? EMPTY_SEQUENCE : bs) :
                    (bs == null ? prefix : new ByteSequence(prefix, bs));
            attributes.add(new KeyValue(key, value));
        } while ( finalByte != '>' );
        needData();
        if ( bufferIterator.next() != '\n' ) {
            throw new MalformedVCFException("unexpected characters at end of metadata line");
        }
        return attributes;
    }

    /** tokenize a tab-delimited data record */
    private List<ByteSequence> captureColumns() {
        final List<ByteSequence> columns = new ArrayList<>();
        byte finalByte;
        do {
            ByteSequence prefix = null;
            int start = bufferIterator.mark();
            do {
                // if we've hit the end of the current buffer, but haven't found a delimiter
                if ( !bufferIterator.hasNext() ) {
                    final ByteSequence bs = bufferIterator.getSequence(start);
                    if ( bs != null ) {
                        prefix = prefix == null ? bs : new ByteSequence(prefix, bs);
                    }
                    needData(); // get a new buffer
                    start = 0;
                }
                finalByte = bufferIterator.next();
            } while ( finalByte != '\t' && finalByte != '\n' );
            final ByteSequence bs = bufferIterator.getSequenceNoDelim(start);
            // if we changed iterators mid-token we'll need to glue two ByteSequences together
            final ByteSequence column = prefix == null ?
                    (bs == null ? EMPTY_SEQUENCE : bs) :
                    (bs == null ? prefix : new ByteSequence(prefix, bs));
            columns.add(column);
        } while ( finalByte != '\n' );
        return columns;
    }

    /** A runtime exception that says that we got confused when reading the vcf:
     * either the format is messed up, or there's an I/O exception. */
    public static final class MalformedVCFException extends RuntimeException {
        public final static long serialVersionUID = 1;

        public MalformedVCFException( final String msg ) {
            super(msg);
        }

        public MalformedVCFException( final String msg, final IOException ioe ) {
            super(msg, ioe);
        }
    }

    /** an iterator over an array of bytes.
     * if it hasNext(), you can do peek() and skip() or you can do next().  it's up to you. */
    public static final class ByteIterator {
        private final byte[] buffer;
        private final int length;
        private int index;

        public ByteIterator( final byte[] buffer, final int start, final int length ) {
            this.buffer = buffer;
            this.length = length;
            this.index = start;
        }

        public boolean hasNext() { return index < length; }

        public byte peek() {
            if ( !hasNext() ) {
                throw new NoSuchElementException();
            }
            return buffer[index];
        }

        public void skip() { index += 1; }

        public byte next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException();
            }
            return buffer[index++];
        }

        /** return the current index.
         * we'll use it as the starting point for a ByteSequence we're constructing */
        public int mark() { return index; }

        /** make a sub-sequence from our buffer from a starting point (that you got by
         * calling mark() earlier) up to the current position */
        public ByteSequence getSequence( final int start ) { return getSequence(start, 0); }

        /** make a sub-sequence ignoring one or more trailing bytes */
        public ByteSequence getSequenceNoDelim( final int start ) { return getSequence(start, 1); }
        private ByteSequence getSequence( final int start, final int skip ) {
            int end = index - skip;
            return end <= start ? null : new ByteSequence(buffer, start, end);
        }
    }

    /** immutable sub-sequence of a byte[] */
    public static final class ByteSequence {
        private final byte[] buffer;
        private final int start;
        private final int end;

        public static final int MISSING_VALUE = Integer.MIN_VALUE;

        public ByteSequence( final byte[] buffer, final int start, final int end ) {
            this.buffer = buffer;
            this.start = start;
            this.end = end;
        }

        public ByteSequence( final String seq ) {
            this.buffer = seq.getBytes();
            this.start = 0;
            this.end = buffer.length;
        }

        public ByteSequence( final ByteSequence seq1, final ByteSequence seq2 ) {
            final int seq1Len = seq1.length();
            final int seq2Len = seq2.length();
            buffer = new byte[seq1Len + seq2Len];
            System.arraycopy(seq1.buffer, seq1.start, buffer, 0, seq1Len);
            System.arraycopy(seq2.buffer, seq2.start, buffer, seq1Len, seq2Len);
            start = 0;
            end = buffer.length;
        }

        public ByteSequence( final ByteSequence... seqs ) {
            int totalLen = 0;
            for ( final ByteSequence seq : seqs ) {
                totalLen += seq.length();
            }
            buffer = new byte[totalLen];
            start = 0;
            end = totalLen;
            int curLen = 0;
            for ( final ByteSequence seq : seqs ) {
                final int len = seq.length();
                System.arraycopy(seq.buffer, seq.start, buffer, curLen, len);
                curLen += len;
            }
        }

        public ByteSequence( final List<ByteSequence> pieces, final char delim ) {
            final int nPieces = pieces.size();
            int totalLen = 0;
            if ( nPieces > 0 ) {
                totalLen = nPieces - 1; // this many delimiters
                for ( final ByteSequence piece : pieces ) {
                    totalLen += piece.length();
                }
            }
            buffer = new byte[totalLen];
            start = 0;
            end = totalLen;
            if ( nPieces > 0 ) {
                ByteSequence piece = pieces.get(0);
                int destIdx = piece.length();
                System.arraycopy(piece.buffer, piece.start, buffer, 0, destIdx);
                for ( int pieceIdx = 1; pieceIdx < nPieces; ++pieceIdx ) {
                    buffer[destIdx++] = (byte)delim;
                    piece = pieces.get(pieceIdx);
                    int len = piece.length();
                    System.arraycopy(piece.buffer, piece.start, buffer, destIdx, len);
                    destIdx += len;
                }
            }
        }

        public int length() { return end - start; }

        public boolean contains( final ByteSequence subSeq ) {
            final int len = subSeq.length();
            final int stop = end - len;
            for ( int idx = start; idx <= stop; ++idx ) {
                int idx1 = idx;
                int idx2 = subSeq.start;
                int nnn = len;
                while ( nnn-- > 0 ) {
                    if ( buffer[idx1++] != subSeq.buffer[idx2++] ) {
                        break;
                    }
                }
                if ( nnn < 0 ) {
                    return true;
                }
            }
            return false;
        }

        public int asInt() {
            final ByteIterator itr = iterator();
            if ( !itr.hasNext() || itr.peek() == '.' ) {
                return MISSING_VALUE;
            }
            int value = 0;
            boolean negative = false;
            if ( itr.peek() == '-' ) {
                negative = true;
                itr.skip();
            }
            while ( itr.hasNext() ) {
                int digit = itr.next();
                if ( digit < '0' || digit > '9' ) {
                    throw new NumberFormatException();
                }
                value = 10 * value + (digit & 0xF);
            }
            return negative ? -value : value;
        }

        public ByteIterator iterator() { return new ByteIterator(buffer, start, end); }

        /** tokenize this ByteSequence on some delimiter */
        public List<ByteSequence> split( final char delim ) {
            final List<ByteSequence> splits = new ArrayList<>();
            final ByteIterator itr = iterator();
            int mark = itr.mark();
            while ( itr.hasNext() ) {
                if ( itr.next() == delim ) {
                    splits.add(itr.getSequenceNoDelim(mark));
                    mark = itr.mark();
                }
            }
            splits.add(itr.getSequence(mark));
            return splits;
        }

        public void write( final OutputStream os ) throws IOException {
            os.write(buffer, start, length());
        }

        @Override public String toString() {
            return new String(buffer, start, length());
        }

        @Override public int hashCode() {
            int result = 101;
            for ( int idx = start; idx < end; ++idx ) {
                result = 47 * result + buffer[idx];
            }
            return 101 * result;
        }

        @Override public boolean equals( final Object obj ) {
            if ( this == obj ) return true;
            if ( !(obj instanceof ByteSequence) ) return false;
            return equals((ByteSequence)obj);
        }

        public boolean equals( final ByteSequence that ) {
            if ( that == null || length() != that.length() ) return false;
            int idx2 = that.start;
            for ( int idx = start; idx < end; ++idx ) {
                if ( buffer[idx] != that.buffer[idx2++] ) return false;
            }
            return true;
        }
    }

    public interface Metadata {
        ByteSequence getKey();
        Object getValue();
        void write( OutputStream os ) throws IOException;
    }

    public static final class KeyValue {
        private final ByteSequence key;
        private final ByteSequence value;

        public KeyValue( final ByteSequence key, final ByteSequence value ) {
            this.key = key;
            this.value = value;
        }

        public ByteSequence getKey() { return key; }
        public ByteSequence getValue() { return value; }

        public void write( final OutputStream os ) throws IOException {
            key.write(os);
            if ( value != null ) {
                os.write('=');
                value.write(os);
            }
        }

        @Override public String toString() { return key + "=" + value; }
    }

    public static final class KeyValueMetadata implements Metadata {
        private final KeyValue keyValue;

        public KeyValueMetadata( final ByteSequence key, final ByteSequence value ) {
            keyValue = new KeyValue(key, value);
        }

        @Override public ByteSequence getKey() { return keyValue.getKey(); }
        @Override public ByteSequence getValue() { return keyValue.getValue(); }

        @Override public void write( final OutputStream os ) throws IOException {
            os.write('#');
            os.write('#');
            keyValue.write(os);
            os.write('\n');
        }

        @Override public String toString() { return keyValue.toString(); }
    }

    public static final class KeyAttributesMetadata implements Metadata {
        private final ByteSequence key;
        private final List<KeyValue> values;

        public KeyAttributesMetadata( final ByteSequence key, final List<KeyValue> values ) {
            this.key = key;
            this.values = values;
        }

        @Override public ByteSequence getKey() { return key; }
        @Override public List<KeyValue> getValue() { return values; }

        @Override public void write( final OutputStream os ) throws IOException {
            os.write('#');
            os.write('#');
            key.write(os);
            os.write('=');
            int prefix = '<';
            for ( final KeyValue kv : values ) {
                os.write(prefix);
                kv.write(os);
                prefix = ',';
            }
            os.write('>');
            os.write('\n');
        }

        @Override public String toString() {
            final StringBuilder sb = new StringBuilder();
            sb.append(key).append("=");
            char prefix = '<';
            for ( final KeyValue kv : values ) {
                sb.append(prefix).append(kv.getKey()).append('=').append(kv.getValue());
                prefix = ',';
            }
            return sb.append(">").toString();
        }
    }

    public static final class ColumnHeaderMetadata implements Metadata {
        private final List<ByteSequence> columns;

        public ColumnHeaderMetadata( final List<ByteSequence> columns ) {
            this.columns = columns;
        }

        @Override public ByteSequence getKey() { return EMPTY_SEQUENCE; }
        @Override public List<ByteSequence> getValue() { return columns; }

        @Override public void write( final OutputStream os ) throws IOException {
            int prefix = '#';
            for ( final ByteSequence col : columns ) {
                os.write(prefix);
                col.write(os);
                prefix = '\t';
            }
            os.write('\n');
        }

        @Override public String toString() {
            final StringBuilder sb = new StringBuilder();
            char prefix = '#';
            for ( final ByteSequence col : columns ) {
                sb.append(prefix).append(col);
                prefix = '\t';
            }
            return sb.toString();
        }
    }

    /** a field like format and genotype with delimited subfields */
    public static final class CompoundField extends AbstractList<ByteSequence> {
        private ByteSequence value;
        private final char delim;
        private List<ByteSequence> subFields;

        public CompoundField( final ByteSequence value, final char delim ) {
            this.value = value;
            this.delim = delim;
            subFields = null;
        }

        public CompoundField( final List<ByteSequence> vals, final char delim ) {
            this.value = null;
            this.delim = delim;
            this.subFields = vals;
        }

        public ByteSequence getValue() {
            if ( value == null ) {
                value = new ByteSequence(subFields, delim);
            }
            return value;
        }

        public void write( final OutputStream os ) throws IOException {
            if ( value != null ) value.write(os);
            else {
                int len = subFields.size();
                if ( len <= 0 ) {
                    os.write('.');
                } else {
                    subFields.get(0).write(os);
                    for ( int idx = 1; idx < len; ++idx ) {
                        os.write(delim);
                        subFields.get(idx).write(os);
                    }
                }
            }
        }

        @Override public int size() {
            populateSubFields();
            return subFields.size();
        }

        @Override public ByteSequence get( final int index ) {
            populateSubFields();
            return subFields.get(index);
        }

        @Override public ByteSequence set( final int index, final ByteSequence val ) {
            populateSubFields();
            value = null;
            return subFields.set(index, val);
        }

        @Override public void add( final int index, final ByteSequence val ) {
            populateSubFields();
            value = null;
            subFields.add(index, val);
        }

        @Override public ByteSequence remove( final int index ) {
            populateSubFields();
            value = null;
            return subFields.remove(index);
        }

        @Override public boolean equals( final Object obj ) {
            if ( this == obj ) return true;
            if ( !(obj instanceof CompoundField) ) return false;
            return getValue().equals(((CompoundField)obj).getValue());
        }
        @Override public int hashCode() {
            return getValue().hashCode();
        }
        @Override public String toString() { return getValue().toString(); }

        private void populateSubFields() {
            if ( subFields == null ) {
                subFields = value.split(delim);
            }
        }
    }

    /** the info subfields are semicolon delimited and contain key/value pairs */
    public static final class InfoField extends AbstractMap<ByteSequence, ByteSequence> {
        private ByteSequence value;
        private LinkedHashMap<ByteSequence, ByteSequence> subFields;

        public InfoField( final ByteSequence value ) {
            this.value = value;
            subFields = null;
        }

        public ByteSequence getValue() {
            if ( value == null ) {
                final ByteArrayOutputStream os = new ByteArrayOutputStream();
                try {
                    write(os);
                } catch ( final IOException ioe ) {
                    throw new IllegalStateException("IOException when writing to ByteArrayOutputStream!?");
                }
                final byte[] buffer = os.toByteArray();
                value = new ByteSequence(buffer, 0, buffer.length);
            }
            return value;
        }

        public void write( final OutputStream os ) throws IOException {
            if ( value != null ) {
                value.write(os);
            } else if ( subFields.isEmpty() ) {
                os.write('.');
            } else {
                boolean needSep = false;
                for ( final Map.Entry<ByteSequence, ByteSequence> entry : subFields.entrySet() ) {
                    if ( needSep ) {
                        os.write(';');
                    }
                    needSep = true;
                    entry.getKey().write(os);
                    final ByteSequence value = entry.getValue();
                    if ( value != null ) {
                        os.write('=');
                        value.write(os);
                    }
                }
            }
        }

        @Override public Set<Entry<ByteSequence, ByteSequence>> entrySet() {
            populateSubFields();
            return subFields.entrySet();
        }

        @Override public boolean containsKey( final Object key ) {
            populateSubFields();
            return subFields.containsKey(key);
        }

        @Override public ByteSequence get( final Object key ) {
            populateSubFields();
            return subFields.get(key);
        }

        @Override public ByteSequence put( final ByteSequence key, final ByteSequence val ) {
            populateSubFields();
            value = null;
            return subFields.put(key, val);
        }

        @Override public ByteSequence remove( final Object key ) {
            populateSubFields();
            if ( containsKey(key) ) {
                value = null;
            }
            return subFields.remove(key);
        }

        private void populateSubFields() {
            if ( subFields == null ) {
                subFields = new LinkedHashMap<>();
                final ByteIterator itr = value.iterator();
                int mark = itr.mark();
                ByteSequence key = null;
                while ( itr.hasNext() ) {
                    byte nextByte = itr.next();
                    if ( nextByte == '=' ) {
                        key = itr.getSequenceNoDelim(mark);
                        mark = itr.mark();
                    } else if ( nextByte == ';' ) {
                        if ( key == null ) {
                            subFields.put(itr.getSequenceNoDelim(mark), null);
                        } else {
                            subFields.put(key, itr.getSequenceNoDelim(mark));
                        }
                        key = null;
                        mark = itr.mark();
                    }
                }
                if ( key == null ) {
                    subFields.put(itr.getSequence(mark), null);
                } else {
                    subFields.put(key, itr.getSequence(mark));
                }
            }
        }
    }

    /** a line of data from the VCF */
    public static final class Record {
        private static final int UNINITIALIZED = -1;

        private final List<ByteSequence> simpleFields;
        private CompoundField filters;
        private InfoField infos;
        private CompoundField formats;
        private final List<CompoundField> genotypes;

        private int position = UNINITIALIZED;
        private int quality = UNINITIALIZED;

        public Record( final List<ByteSequence> vals ) {
            simpleFields = new ArrayList<>(vals.subList(0, 6));
            filters = new CompoundField(vals.get(6), ';');
            infos = new InfoField(vals.get(7));
            final int nVals = vals.size();
            formats = nVals > 8 ? new CompoundField(vals.get(8), ':') : null;
            genotypes = new ArrayList<>(Math.max(0, nVals - 9));
            for ( int idx = 9; idx < nVals; ++idx ) {
                genotypes.add(new CompoundField(vals.get(idx), ':'));
            }
        }

        public ByteSequence getChromosome() { return simpleFields.get(0); }
        public void setChromosome( final ByteSequence val ) { simpleFields.set(0, val); }

        public int getPosition() {
            if ( position == UNINITIALIZED ) {
                position = simpleFields.get(1).asInt();
            }
            return position;
        }
        public void setPosition( final int pos ) {
            setPosition(new ByteSequence(Integer.toString(pos)));
        }
        public void setPosition( final ByteSequence val ) {
            simpleFields.set(1, val);
            position = UNINITIALIZED;
        }

        public ByteSequence getID() { return simpleFields.get(2); }
        public void setID( final ByteSequence val ) { simpleFields.set(2, val); }

        public ByteSequence getRef() { return simpleFields.get(3); }
        public void setRef( final ByteSequence val ) { simpleFields.set(3, val); }

        public ByteSequence getAlt() { return simpleFields.get(4); }
        public void setAlt( final ByteSequence val ) { simpleFields.set(4, val); }

        public int getQuality() {
            if ( quality == UNINITIALIZED ) {
                quality = simpleFields.get(5).asInt();
            }
            return quality;
        }
        public void setQuality( final ByteSequence val ) {
            simpleFields.set(5, val);
            quality = UNINITIALIZED;
        }

        public CompoundField getFilter() { return filters; }
        public void setFilter( final ByteSequence val ) {
            filters = new CompoundField(val, ';');
        }
        public void setFilter( final List<ByteSequence> vals ) {
            filters = new CompoundField(vals, ';');
        }

        public InfoField getInfo() {
            return infos;
        }
        public void setInfo( final ByteSequence val ) { infos = new InfoField(val); }

        public CompoundField getFormat() { return formats; }
        public void setFormat( final ByteSequence val ) {
            formats = new CompoundField(val, ':');
        }

        public List<CompoundField> getGenotypes() { return genotypes; }
        public void setGenotypes( final List<ByteSequence> vals ) {
            genotypes.clear();
            for ( final ByteSequence val : vals ) {
                genotypes.add(new CompoundField(val, ':'));
            }
        }

        public void write( final OutputStream os ) throws IOException {
            simpleFields.get(0).write(os);
            for ( int idx = 1; idx < 6; ++idx ) {
                os.write('\t');
                simpleFields.get(idx).write(os);
            }
            os.write('\t');
            filters.write(os);
            os.write('\t');
            infos.write(os);
            if ( formats != null ) {
                os.write('\t');
                formats.write(os);
                for ( final CompoundField genotype : genotypes ) {
                    os.write('\t');
                    genotype.write(os);
                }
            }
            os.write('\n');
        }

        @Override
        public String toString() {
            final StringBuilder sb = new StringBuilder();
            String prefix = "";
            for ( final ByteSequence field : simpleFields ) {
                sb.append(prefix).append(field.toString());
                prefix = "\t";
            }
            return sb.toString();
        }
    }
}
