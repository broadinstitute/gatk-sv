package org.broadinstitute.svpipeline;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.svpipeline.VCFParser.*;

public final class VCFParserUnitTest {
    public static void main( final String[] args ) {
        testAsserts();
        testEmptyFile();
        testFileFormatMetadata();
        testFilter();
        testColumnHeaders();
        testRecord();
        testRoundTrip();
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
            throw new AssertionError("assertions aren't turned on (with -ea), so you're not testing anything.");
        }
    }

    public static void testEmptyFile() {
        boolean caughtIt = false;
        try ( final VCFParser parser = new VCFParser("/dev/null") ) {
            assert(!parser.hasMetadata());
        } catch ( final MalformedVCFException emptyVCF ) {
            caughtIt = true;
        }
        assert(caughtIt);
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
        parser.close();
    }

    public static void testFilter() {
        final byte[] bytes = "##FILTER=<ID=PASS,Description=\"All filters passed\">\n".getBytes();
        final VCFParser parser = new VCFParser(new ByteArrayInputStream(bytes));
        assert(parser.hasMetadata());
        final Metadata metadata = parser.nextMetaData();
        assert(metadata instanceof KeyAttributesMetadata);
        final KeyAttributesMetadata kaMetadata = (KeyAttributesMetadata)metadata;
        assert(kaMetadata.getKey().equals(new ByteSequence("FILTER")));
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
        final String line = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2";
        final byte[] bytes = ("#" + line + "\n").getBytes();
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

    public static void testRecord() {
        final String line = "chr1\t10000\tna19240_DUP_chr1_1\tN\t<DUP>\t999\tPASS;BUT_FUNKY\t" +
                "END=16000;SVTYPE=DUP;FLAG1;CHR2=chr1;SVLEN=6000;ALGORITHMS=depth;EVIDENCE=RD;FLAG2\t" +
                "GT:GQ:RD_CN:RD_GQ:PE_GT:PE_GQ:SR_GT:SR_GQ:EV\t0/1:142:3:142:.:.:.:.:RD\t" +
                "0/0:999:2:999:.:.:.:.:RD";
        final byte[] bytes = (line + "\n").getBytes();
        final VCFParser parser = new VCFParser(new ByteArrayInputStream(bytes));
        assert(!parser.hasMetadata());
        assert(parser.hasRecord());
        final Record record = parser.nextRecord();
        assert(!parser.hasMetadata());
        assert(!parser.hasRecord());
        final String[] cols = line.split("\t");

        assert(record.getChromosome().equals(new ByteSequence(cols[0])));
        final ByteSequence newChr = new ByteSequence("chr1");
        record.setChromosome(newChr);
        assert(record.getChromosome().equals(newChr));

        final int curPosition = record.getPosition();
        assert(curPosition == Integer.parseInt(cols[1]));
        final ByteSequence newPos = new ByteSequence("10001");
        record.setPosition(newPos);
        final int newPosition = record.getPosition();
        assert(newPosition == 10001);

        assert(record.getID().equals(new ByteSequence(cols[2])));
        final ByteSequence newID = new ByteSequence("newID");
        record.setID(newID);
        assert(record.getID().equals(newID));

        assert(record.getRef().equals(new ByteSequence(cols[3])));
        final ByteSequence newRef = new ByteSequence("A");
        record.setRef(newRef);
        assert(record.getRef().equals(newRef));

        assert(record.getAlt().equals(new ByteSequence(cols[4])));
        final ByteSequence newAlt = new ByteSequence("C");
        record.setAlt(newAlt);
        assert(record.getAlt().equals(newAlt));

        final int curQuality = record.getQuality();
        assert(curQuality == Integer.parseInt(cols[5]));
        final ByteSequence newQual = new ByteSequence("1");
        record.setQuality(newQual);
        final int newQuality = record.getQuality();
        assert(newQuality == 1);

        final CompoundField filters = record.getFilter();
        final ByteSequence originalFilters = new ByteSequence(cols[6]);
        final ByteSequence curFilters = filters.getValue();
        assert(curFilters.equals(originalFilters));
        assert(filters.size() == 2);
        assert(filters.get(0).equals(new ByteSequence("PASS")));
        assert(filters.get(1).equals(new ByteSequence("BUT_FUNKY")));
        final ByteSequence failFilter = new ByteSequence("FAIL");
        filters.set(0, failFilter);
        assert(filters.get(0).equals(failFilter));
        final ByteSequence newFilters = filters.getValue();
        assert(newFilters.equals(new ByteSequence("FAIL;BUT_FUNKY")));
        record.setFilter(originalFilters);
        final CompoundField revisedFilters = record.getFilter();
        final ByteSequence newerFilters = revisedFilters.getValue();
        assert(newerFilters.equals(originalFilters));
        revisedFilters.add(revisedFilters.remove(0));
        final ByteSequence newestFilters = revisedFilters.getValue();
        assert(newestFilters.equals(new ByteSequence("BUT_FUNKY;PASS")));

        final InfoField info = record.getInfo();
        final ByteSequence originalInfo = new ByteSequence(cols[7]);
        final ByteSequence curInfo = info.getValue();
        assert(curInfo.equals(originalInfo));
        final String[] infoVals = cols[7].split(";");
        assert(info.size() == infoVals.length);
        for ( final String val : infoVals ) {
            final String[] kv = val.split("=");
            final ByteSequence key = new ByteSequence(kv[0]);
            assert(info.containsKey(key));
            if ( kv.length > 1 ) {
                assert(info.get(key).equals(new ByteSequence(kv[1])));
            } else {
                assert(info.get(key) == null);
            }
        }
        final ByteSequence svLenKey = new ByteSequence("SVLEN");
        final ByteSequence newSVLen = new ByteSequence("6001");
        info.put(svLenKey, newSVLen);
        assert(info.get(svLenKey).equals(newSVLen));
        info.put(svLenKey, new ByteSequence("6000"));
        final ByteSequence newInfoValue = info.getValue();
        assert(newInfoValue.equals(originalInfo));
        final ByteSequence flag1Key = new ByteSequence("FLAG1");
        info.remove(flag1Key);
        assert(info.get(flag1Key) == null);
        final ByteSequence flag2Key = new ByteSequence("FLAG2");
        record.setInfo(flag2Key);
        final InfoField newInfo = record.getInfo();
        assert(!newInfo.containsKey(flag1Key));
        assert(newInfo.containsKey(flag2Key));

        final CompoundField format = record.getFormat();
        final ByteSequence originalFormat = new ByteSequence(cols[8]);
        final ByteSequence curFormat = format.getValue();
        assert(curFormat.equals(originalFormat));
        record.setFormat(new ByteSequence("GT"));
        assert(record.getFormat().size() == 1);

        final List<CompoundField> genotypes = record.getGenotypes();
        assert(genotypes.size() == 2);
        final ByteSequence geno1 = genotypes.get(0).getValue();
        final ByteSequence geno1Value = new ByteSequence("0/1:142:3:142:.:.:.:.:RD");
        assert(geno1.equals(geno1Value));
        final ByteSequence geno2 = genotypes.get(1).getValue();
        final ByteSequence geno2Value = new ByteSequence("0/0:999:2:999:.:.:.:.:RD");
        assert(geno2.equals(geno2Value));
        record.setGenotypes(Collections.singletonList(geno2Value));
        final List<CompoundField> newGenotypes = record.getGenotypes();
        assert(newGenotypes.size() == 1);
        final ByteSequence newGeno2Value = newGenotypes.get(0).getValue();
        assert(newGeno2Value.equals(geno2Value));
    }

    public static void testRoundTrip() {
        final StringBuilder sb = new StringBuilder(100000);
        for ( int idx = 0; idx < 1000; ++idx ) {
            buildLine(idx, sb);
        }
        final byte[] bytes = sb.toString().getBytes();
        final ByteArrayInputStream is = new ByteArrayInputStream(bytes);
        final VCFParser parser = new VCFParser(is);
        final ByteArrayOutputStream os = new ByteArrayOutputStream(100000);
        while ( parser.hasRecord() ) {
            final Record record = parser.nextRecord();
            try {
                record.write(os);
            } catch ( final IOException ioe ) {
                throw new RuntimeException("unexpected IOException");
            }
        }
        parser.close();
        assert(Arrays.equals(os.toByteArray(),bytes));
    }

    private static void buildLine( final int idx, final StringBuilder sb ) {
        final int pos = 10000 * idx;
        sb.append("chr1\t").append(10000+100*idx).append('\t').append("Event").append(idx).append('\t');
        sb.append("N\t").append("<DUP>\t").append("999\t").append("PASS\t");
        sb.append("END=").append(pos+999).append('\t');
        sb.append("SVTYPE=DUP;CHR2=chr1;SVLEN=1000;ALGORITHMS=depth;EVIDENCE=RD\t");
        sb.append("GT:GQ:RD_CN\t").append("0/1:999:2\t").append("0/0:999:1\n");
    }
}
