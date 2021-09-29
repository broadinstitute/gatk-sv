package org.broadinstitute.svpipeline;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.broadinstitute.svpipeline.VCFParser.*;

public final class VCFParserUnitTest {
    public static void main( final String[] args ) {
        testFileFormatMetadata();
        testFilter();
        testColumnHeaders();
        System.out.println("OK");
    }

    public static void testFileFormatMetadata() {
        final byte[] bytes = "##fileformat=VCFv4.2\n".getBytes();
        final VCFParser parser = new VCFParser(new ByteArrayInputStream(bytes));
        assert(parser.hasMetadata());
        final Metadata metadata = parser.nextMetaData();
        assert(metadata instanceof KeyValueMetadata);
        final KeyValueMetadata kvMetadata = (KeyValueMetadata)metadata;
        assert(kvMetadata.getKey().equals(new ByteSequence("fileformat")));
        assert(kvMetadata.getValue().equals(new ByteSequence("VCFv4.2")));
        assert(!parser.hasMetadata());
        assert(!parser.hasRecord());
        try ( final ByteArrayOutputStream os = new ByteArrayOutputStream() ) {
            metadata.write(os);
            assert(Arrays.equals(bytes, os.toByteArray()));
        } catch ( final IOException ioe ) {
            throw new RuntimeException(ioe);
        }
    }

    public static void testFilter() {
        final byte[] bytes = "##FILTER=<ID=PASS,Description=\"All filters passed\">\n".getBytes();
        final VCFParser parser = new VCFParser(new ByteArrayInputStream(bytes));
        assert(parser.hasMetadata());
        final Metadata metadata = parser.nextMetaData();
        assert(metadata instanceof KeyAttributesMetadata);
        final KeyAttributesMetadata kaMetadata = (KeyAttributesMetadata)metadata;
        assert(kaMetadata.getKey().equals(new ByteSequence("FORMAT")));
        final List<KeyValue> kaValues = kaMetadata.getValue();
        assert(kaValues.size() == 2);
        final KeyValue kv0 = kaValues.get(0);
        assert(kv0.getKey().equals(new ByteSequence("ID")));
        assert(kv0.getValue().equals(new ByteSequence("PASS")));
        final KeyValue kv1 = kaValues.get(1);
        assert(kv1.getKey().equals(new ByteSequence("Description")));
        assert(kv1.getValue().equals(new ByteSequence("\"All filters passed\"")));
        assert(!parser.hasMetadata());
        assert(!parser.hasRecord());
        try ( final ByteArrayOutputStream os = new ByteArrayOutputStream() ) {
            metadata.write(os);
            assert(Arrays.equals(bytes, os.toByteArray()));
        } catch ( final IOException ioe ) {
            throw new RuntimeException(ioe);
        }
    }

    public static void testColumnHeaders() {
        final String line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\n";
        final byte[] bytes = line.getBytes();
        final VCFParser parser = new VCFParser(new ByteArrayInputStream(bytes));
        assert(parser.hasMetadata());
        final Metadata metadata = parser.nextMetaData();
        assert(metadata instanceof ColumnHeaderMetadata);
        final ColumnHeaderMetadata columns = (ColumnHeaderMetadata)metadata;
        final List<ByteSequence> cols = columns.getValue();
        final String[] splitLine = line.split("\t");
        assert(splitLine.length == cols.size());
        for ( int idx = 0; idx < splitLine.length; ++idx ) {
            assert(cols.get(idx).equals(new ByteSequence(splitLine[idx])));
        }
        assert(!parser.hasMetadata());
        assert(!parser.hasRecord());
        try ( final ByteArrayOutputStream os = new ByteArrayOutputStream() ) {
            metadata.write(os);
            assert(Arrays.equals(bytes, os.toByteArray()));
        } catch ( final IOException ioe ) {
            throw new RuntimeException(ioe);
        }
    }
}
