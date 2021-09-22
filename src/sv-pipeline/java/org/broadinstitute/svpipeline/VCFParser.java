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
            this.is = System.in instanceof BufferedInputStream ?
                    System.in :
                    new BufferedInputStream(System.in);
        } else {
            this.pathName = pathName;
            try {
                final BufferedInputStream bis =
                        new BufferedInputStream(new FileInputStream(pathName));
                this.is = pathName.endsWith(GZ) ? new GZIPInputStream(bis) : bis;
            } catch ( final IOException ioe ) {
                throw new MalformedVCFException("can't open " + pathName, ioe);
            }
        }
        if ( !readBuffer() ) {
            throw new MalformedVCFException(this.pathName + " is empty");
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
            return new Columns(captureColumns());
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
            return new KeyValue(key, value);
        }
        bufferIterator.skip();
        // yup.  multiple values. tokenize them.
        return new KeyAttributes(key, captureAttributes());
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

    private void expect( final String expect ) {
        final int expectLen = expect.length();
        for ( int iii = 0; iii < expectLen; ++iii ) {
            needData();
            final byte nextByte = bufferIterator.next();
            if ( expect.charAt(iii) != nextByte ) {
                throw new MalformedVCFException("expected " + expect + " but found " +
                        expect.substring(0, iii) + (char)nextByte + "...");
            }
        }
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
        expect("\n");
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
        public static final ByteSequence EMPTY = new ByteSequence(new byte[0], 0, 0);

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

        public ByteSequence replace( final ByteSequence oldValue, final ByteSequence newValue ) {
            if ( buffer != oldValue.buffer || oldValue.start < start || oldValue.end > end ) {
                throw new IllegalStateException("replaced value not within buffer");
            }
            final int oldLen = oldValue.length();
            final int newLen = newValue.length();
            if ( newLen == oldLen ) {
                System.arraycopy(newValue.buffer, newValue.start, buffer, oldValue.start, newLen);
                System.arraycopy(buffer, oldValue.end, buffer, oldValue.start + newLen, end - oldValue.end);
                return this;
            }
            final int length = length();
            final byte[] newBuf = new byte[length + newLen - oldValue.length()];
            final int len1 = oldValue.start - start;
            System.arraycopy(buffer, start, newBuf, 0, len1);
            System.arraycopy(newValue.buffer, newValue.start, newBuf, len1, newLen);
            System.arraycopy(buffer, oldValue.end, newBuf, len1 + newLen, end - oldValue.end);
            return new ByteSequence(newBuf, 0, newBuf.length);
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

    enum MetadataType {
        KeyValue,
        KeyAttributes,
        Columns
    }

    public interface Metadata {
        MetadataType getType();
        ByteSequence getKey();
        Object getValue();
        void write( OutputStream os ) throws IOException;
    }

    public static final class KeyValue implements Metadata {
        private final ByteSequence key;
        private final ByteSequence value;

        public KeyValue( final ByteSequence key, final ByteSequence value ) {
            this.key = key;
            this.value = value;
        }

        @Override public MetadataType getType() { return MetadataType.KeyValue; }
        @Override public ByteSequence getKey() { return key; }
        @Override public ByteSequence getValue() { return value; }

        @Override public void write( final OutputStream os ) throws IOException {
            os.write('#');
            os.write('#');
            key.write(os);
            if ( value != null ) {
                os.write('=');
                value.write(os);
            }
            os.write('\n');
        }

        @Override public String toString() { return key + "=" + value; }
    }

    public static final class KeyAttributes implements Metadata {
        private final ByteSequence key;
        private final List<KeyValue> values;

        public KeyAttributes( final ByteSequence key, final List<KeyValue> values ) {
            this.key = key;
            this.values = values;
        }

        @Override public MetadataType getType() { return MetadataType.KeyAttributes; }
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
                kv.getKey().write(os);
                os.write('=');
                kv.getValue().write(os);
                prefix = ',';
            }
            os.write('>');
            os.write('\n');
        }

        @Override public String toString() {
            final StringBuilder sb = new StringBuilder();
            sb.append("##").append(key).append("=");
            char prefix = '<';
            for ( final KeyValue kv : values ) {
                sb.append(prefix).append(kv.getKey()).append('=').append(kv.getValue());
                prefix = ',';
            }
            return sb.append(">").toString();
        }
    }

    public static final class Columns implements Metadata {
        private final List<ByteSequence> columns;

        public Columns( final List<ByteSequence> columns ) {
            this.columns = columns;
        }

        @Override public MetadataType getType() { return MetadataType.Columns; }
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

    /** a line of data from the VCF */
    public static final class Record {
        private static final int UNINITIALIZED = -1;

        private final List<ByteSequence> columns;
        private List<KeyValue> infoKeyValues = null;
        private int position = UNINITIALIZED;
        private int quality = UNINITIALIZED;

        public Record( final List<ByteSequence> columns ) {
            this.columns = columns;
        }

        public ByteSequence getChromosome() { return columns.get(0); }

        public int getPosition() {
            if ( position == UNINITIALIZED ) {
                position = columns.get(1).asInt();
            }
            return position;
        }

        public ByteSequence getID() { return columns.get(2); }
        public ByteSequence getRef() { return columns.get(3); }
        public void setRef( final ByteSequence val ) { columns.set(3, val); }
        public ByteSequence getAlt() { return columns.get(4); }
        public void setAlt( final ByteSequence val ) { columns.set(4, val); }

        public int getQuality() {
            if ( quality == UNINITIALIZED ) {
                quality = columns.get(5).asInt();
            }
            return quality;
        }
        public void setQuality( final ByteSequence val ) { columns.set(5, val); }

        public ByteSequence getFilter() { return columns.get(6); }
        public void setFilter( final ByteSequence val ) { columns.set(6, val); }

        public List<KeyValue> getInfo() {
            if ( infoKeyValues == null ) {
                infoKeyValues = parseKVs(columns.get(7));
            }
            return infoKeyValues;
        }

        public ByteSequence getInfoField( final ByteSequence key ) {
            final List<KeyValue> kvs = getInfo();
            for ( final KeyValue kv : kvs ) {
                if ( key.equals(kv.getKey()) ) return kv.getValue();
            }
            return null;
        }

        public void setInfoField( final ByteSequence oldValue, final ByteSequence newValue ) {
            infoKeyValues = null;
            columns.set(7, columns.get(7).replace(oldValue, newValue));
        }

        public void removeInfoField( final ByteSequence key ) {
            final List<KeyValue> kvs = getInfo();
            for ( final KeyValue kv : kvs ) {
                final ByteSequence kvKey = kv.getKey();
                if ( key.equals(kvKey) ) {
                    ByteSequence info = columns.get(7);
                    int start = kvKey.start;
                    final ByteSequence kvValue = kv.getValue();
                    int end = kvValue == null ? kvKey.end : kvValue.end;
                    if ( start > info.start ) start -= 1;
                    else if ( end < info.end ) end += 1;
                    columns.set(7, info.replace(new ByteSequence(info.buffer, start, end), ByteSequence.EMPTY));
                    infoKeyValues = null;
                    break;
                }
            }
        }

        public Map<ByteSequence, ByteSequence> getInfoAsMap() {
            final List<KeyValue> infoList = getInfo();
            final Map<ByteSequence, ByteSequence> infoMap = new HashMap<>(infoList.size() * 2);
            infoList.forEach(kv -> infoMap.put(kv.getKey(), kv.getValue()));
            return infoMap;
        }

        public ByteSequence getFormat() { return columns.size() > 8 ? columns.get(8) : null; }
        public int getFormatIndex( final ByteSequence val ) {
            if ( columns.size() <= 8 ) return -1;
            final List<ByteSequence> formats = columns.get(8).split(':');
            final int nFormats = formats.size();
            for ( int fmtIdx = 0; fmtIdx < nFormats; ++fmtIdx ) {
                final ByteSequence fmt = formats.get(fmtIdx);
                if ( fmt.equals(val) ) {
                    return fmtIdx;
                }
            }
            return -1;
        }

        public List<ByteSequence> getGenotypes() {
            return columns.size() > 9 ? columns.subList(9, columns.size()) : Collections.emptyList();
        }

        public void write( final OutputStream os ) throws IOException {
            final int nCols = columns.size();
            columns.get(0).write(os);
            for ( int iii = 1; iii < nCols; ++iii ) {
                os.write('\t');
                columns.get(iii).write(os);
            }
            os.write('\n');
        }

        @Override public String toString() {
            final StringBuilder sb = new StringBuilder();
            String prefix = "";
            for ( final ByteSequence col : columns ) {
                sb.append(prefix).append(col);
                prefix = "\t";
            }
            return sb.toString();
        }

        private static List<KeyValue> parseKVs( final ByteSequence bs ) {
            final List<KeyValue> attributes = new ArrayList<>();
            final ByteIterator itr = bs.iterator();
            int mark = itr.mark();
            ByteSequence key = null;
            while ( itr.hasNext() ) {
                byte nextByte = itr.next();
                if ( nextByte == '=' ) {
                    key = itr.getSequenceNoDelim(mark);
                    mark = itr.mark();
                } else if ( nextByte == ';' ) {
                    if ( key == null ) {
                        attributes.add(new KeyValue(itr.getSequenceNoDelim(mark), null));
                    } else {
                        attributes.add(new KeyValue(key, itr.getSequenceNoDelim(mark)));
                    }
                    key = null;
                    mark = itr.mark();
                }
            }
            if ( key == null ) {
                attributes.add(new KeyValue(itr.getSequence(mark), null));
            } else {
                attributes.add(new KeyValue(key, itr.getSequence(mark)));
            }
            return attributes;
        }
    }
}
