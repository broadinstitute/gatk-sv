package org.broadinstitute.svpipeline;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.Pattern;
import org.broadinstitute.svpipeline.VCFParser.*;

public class CleanVCFPart1 {
    private static final ByteSequence[] EV_VALS = {
            null,
            new ByteSequence("RD"),
            new ByteSequence("PE"),
            new ByteSequence("RD,PE"),
            new ByteSequence("SR"),
            new ByteSequence("RD,SR"),
            new ByteSequence("PE,SR"),
            new ByteSequence("RD,PE,SR")
    };
    private static final ByteSequence FORMAT_LINE = new ByteSequence("FORMAT");
    private static final ByteSequence ID_KEY = new ByteSequence("ID");
    private static final ByteSequence EV_VALUE = new ByteSequence("EV");
    private static final ByteSequence TYPE_KEY = new ByteSequence("Type");
    private static final ByteSequence STRING_VALUE = new ByteSequence("String");
    private static final ByteSequence NUMBER_KEY = new ByteSequence("Number");
    private static final ByteSequence SVTYPE_KEY = new ByteSequence("SVTYPE");
    private static final ByteSequence ME_VALUE = new ByteSequence(":ME");
    private static final ByteSequence LT_VALUE = new ByteSequence("<");
    private static final ByteSequence GT_VALUE = new ByteSequence(">");
    private static final ByteSequence N_VALUE = new ByteSequence("N");
    private static final ByteSequence END_KEY = new ByteSequence("END");
    private static final ByteSequence VARGQ_KEY = new ByteSequence("varGQ");
    private static final ByteSequence MULTIALLELIC_KEY = new ByteSequence("MULTIALLELIC");
    private static final ByteSequence UNRESOLVED_KEY = new ByteSequence("UNRESOLVED");
    private static final ByteSequence HIGH_SR_BACKGROUND = new ByteSequence("HIGH_SR_BACKGROUND");
    private static final ByteSequence PASS_VALUE = new ByteSequence("PASS");
    private static final ByteSequence BOTHSIDES_VALUE = new ByteSequence("BOTHSIDES_SUPPORT");
    private static final ByteSequence DEL_VALUE = new ByteSequence("DEL");
    private static final ByteSequence DUP_VALUE = new ByteSequence("DUP");
    private static final ByteSequence RDCN_VALUE = new ByteSequence("RD_CN");
    private static final ByteSequence MISSING_VALUE = new ByteSequence(".");
    private static final ByteSequence MISSING_GENOTYPE = new ByteSequence("./.");
    private static final ByteSequence GT_REF_REF = new ByteSequence("0/0");
    private static final ByteSequence GT_REF_ALT = new ByteSequence("0/1");
    private static final ByteSequence GT_ALT_ALT = new ByteSequence("1/1");

    private static final int MIN_ALLOSOME_EVENT_SIZE = 5000;

    public static void main( final String[] args ) {
        if ( args.length != 8 ) {
            System.err.println("Usage: java org.broadinstitute.svpipeline.CleanVCFPart1 " +
                    "INPUTVCFFILE PEDIGREES XCHR YCHR NOISYEVENTS BOTHSIDES SAMPLESOUT REVISEDEVENTSOUT");
            System.exit(1);
        }
        final VCFParser parser = new VCFParser(args[0]);
        final ByteSequence xChrName = new ByteSequence(args[2]);
        final ByteSequence yChrName = new ByteSequence(args[3]);
        final Set<ByteSequence> noisyEvents = readNoisyEventsFile(args[4]);
        final Set<ByteSequence> bothsidesSupportEvents = readBothSidesFile(args[5]);
        try ( final OutputStream os
                      = new BufferedOutputStream(new FileOutputStream(FileDescriptor.out));
              final OutputStream osSamples = new BufferedOutputStream(new FileOutputStream(args[6]));
              final OutputStream osRevEvents = new BufferedOutputStream(new FileOutputStream(args[7])) ) {
            int[] sexForSample = null;
            while ( parser.hasMetadata() ) {
                final Metadata metadata = parser.nextMetaData();
                if ( metadata instanceof ColumnHeaderMetadata ) {
                    final ColumnHeaderMetadata cols = ((ColumnHeaderMetadata)metadata);
                    final List<ByteSequence> colNames = cols.getValue();
                    final int nCols = colNames.size();
                    for ( int idx = 9; idx < nCols; ++idx ) {
                        colNames.get(idx).write(osSamples);
                        osSamples.write('\n');
                    }
                    sexForSample = readPedFile(args[1], cols.getValue());
                    os.write(("##FILTER=<ID=HIGH_SR_BACKGROUND,Description=\"High number of "
                            + "SR splits in background samples indicating messy region\">\n")
                                .getBytes(StandardCharsets.UTF_8));
                    os.write("##FILTER=<ID=UNRESOLVED,Description=\"Variant is unresolved\">\n"
                                .getBytes(StandardCharsets.UTF_8));
                    os.write(("##FILTER=<ID=BOTHSIDES_SUPPORT,Description=\"Variant has " +
                            "read-level support for both sides of breakpoint\">\n")
                                .getBytes(StandardCharsets.UTF_8));
                } else if ( metadata instanceof KeyAttributesMetadata ) {
                    final KeyAttributesMetadata keyAttrs = (KeyAttributesMetadata)metadata;
                    if ( keyAttrs.getKey().equals(FORMAT_LINE) ) {
                        final List<KeyValue> kvs = keyAttrs.getValue();
                        final int nKVs = kvs.size();
                        if ( nKVs > 2 ) {
                            final KeyValue kv0 = kvs.get(0);
                            final KeyValue kv1 = kvs.get(1);
                            final KeyValue kv2 = kvs.get(2);
                            if ( kv0.getKey().equals(ID_KEY) && kv0.getValue().equals(EV_VALUE) ) {
                                if ( kv1.getKey().equals(NUMBER_KEY) ) {
                                    kvs.set(1, new KeyValue(NUMBER_KEY, MISSING_VALUE));
                                }
                                if ( kv2.getKey().equals(TYPE_KEY) ) {
                                    kvs.set(2, new KeyValue(TYPE_KEY, STRING_VALUE));
                                }
                            }
                        }
                    }
                }
                metadata.write(os);
            }
            if ( sexForSample == null ) {
                throw new RuntimeException("header line with sample names is missing.");
            }
            while ( parser.hasRecord() ) {
                final Record record = parser.nextRecord();

                // replace the numeric EV value with a text value
                final int evIdx = record.getFormat().indexOf(EV_VALUE);
                if ( evIdx >= 0 ) {
                    for ( final CompoundField genotypeVals : record.getGenotypes() ) {
                        genotypeVals.set(evIdx, EV_VALS[genotypeVals.get(evIdx).asInt()]);
                    }
                }

                // move the SVTYPE to the ALT field (except for MEs)
                final InfoField info = record.getInfo();
                final ByteSequence svType = info.get(SVTYPE_KEY);
                if ( !record.getAlt().contains(ME_VALUE) ) {
                    if ( svType != null ) {
                        record.setAlt(new ByteSequence(LT_VALUE, svType, GT_VALUE));
                    }
                }
                record.setRef(N_VALUE);

                // move varGQ info field to quality column
                final ByteSequence varGQ = info.get(VARGQ_KEY);
                if ( varGQ != null ) {
                    record.setQuality(varGQ);
                    info.remove(VARGQ_KEY);
                }

                // remove MULTIALLELIC flag, if present
                info.remove(MULTIALLELIC_KEY);

                // remove UNRESOLVED flag and add it as a filter
                if ( info.containsKey(UNRESOLVED_KEY) ) {
                    record.getFilter().add(UNRESOLVED_KEY);
                    info.remove(UNRESOLVED_KEY);
                }

                // mark noisy events
                if ( noisyEvents.contains(record.getID()) ) {
                    record.getFilter().add(HIGH_SR_BACKGROUND);
                }

                // mark bothsides support
                if ( bothsidesSupportEvents.contains(record.getID()) ) {
                    final CompoundField filters = record.getFilter();
                    if ( filters.size() == 1 && filters.get(0).equals(PASS_VALUE) ) {
                        record.setFilter(BOTHSIDES_VALUE);
                    } else {
                        filters.add(BOTHSIDES_VALUE);
                    }
                }

                // fix genotypes on allosomes
                final boolean isY;
                if ( (isY = yChrName.equals(record.getChromosome())) ||
                        xChrName.equals(record.getChromosome())) {
                    final List<CompoundField> genotypes = record.getGenotypes();
                    final int rdCNIndex = record.getFormat().indexOf(RDCN_VALUE);
                    final ByteSequence end = info.get(END_KEY);
                    boolean adjustMale = false;
                    final boolean isDel;
                    if ( ((isDel = DEL_VALUE.equals(svType)) || DUP_VALUE.equals(svType)) && rdCNIndex >= 0 && end != null &&
                            end.asInt() + 1 - record.getPosition() > MIN_ALLOSOME_EVENT_SIZE ) {
                        adjustMale = isRevisableEvent(genotypes, rdCNIndex, sexForSample, isY);
                        if ( adjustMale ) {
                            record.getID().write(osRevEvents);
                            osRevEvents.write('\n');
                        }
                    }
                    CompoundField emptyGenotype = null;
                    final int nSamples = genotypes.size();
                    for ( int sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx ) {
                        final int sampleSex = sexForSample[sampleIdx];
                        final CompoundField genotype = genotypes.get(sampleIdx);
                        if ( sampleSex == 1 ) {
                            if ( adjustMale ) {
                                final ByteSequence rdCN = genotype.get(rdCNIndex);
                                if ( rdCN.equals(MISSING_VALUE) ) {
                                    continue;
                                }
                                final int rdCNVal = rdCN.asInt();
                                genotype.set(rdCNIndex, new ByteSequence(Integer.toString(rdCNVal + 1)));
                                if ( isDel ) {
                                    if ( rdCNVal >= 1 ) genotype.set(0, GT_REF_REF);
                                    else if ( rdCNVal == 0 ) genotype.set(0, GT_REF_ALT);
                                } else {
                                    if ( rdCNVal <= 1 ) genotype.set(0, GT_REF_REF);
                                    else if ( rdCNVal == 2 ) genotype.set(0, GT_REF_ALT);
                                    else genotype.set(0, GT_ALT_ALT);
                                }
                            }
                        } else if ( sampleSex == 2 ) {
                            if ( isY ) {
                                if ( emptyGenotype == null ) {
                                    emptyGenotype = new CompoundField(MISSING_GENOTYPE, ':');
                                    int nFields = genotype.size();
                                    while ( --nFields > 0 ) {
                                        emptyGenotype.add(MISSING_VALUE);
                                    }
                                    emptyGenotype.getValue(); // performance hack to put the pieces together
                                }
                                genotypes.set(sampleIdx, emptyGenotype);
                            }
                        } else {
                            genotype.set(0, MISSING_GENOTYPE);
                        }
                    }
                }

                record.write(os);
            }
        } catch ( final IOException ioe ) {
            throw new RuntimeException("Can't write to stdout", ioe);
        }
    }

    private static boolean isRevisableEvent( final List<CompoundField> genotypes,
                                             final int rdCNIndex,
                                             final int[] sexForColumn,
                                             final boolean isY ) {
        // We're going to calculate the median rdCN values for males and females.
        // We only care if the median is 0, 1, 2, or something larger, so we'll use 4 bins to
        // sum up the counts:  all values >2 go into the last bucket.
        final int[] maleCounts = new int[4];
        final int[] femaleCounts = new int[4];
        final int nSamples = genotypes.size();
        for ( int sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx ) {
            final ByteSequence rdCN = genotypes.get(sampleIdx).get(rdCNIndex);
            if ( MISSING_VALUE.equals(rdCN) ) {
                continue;
            }
            int rdCNVal = rdCN.asInt();
            if ( rdCNVal > 2 ) {
                rdCNVal = 3;
            }
            final int sampleSex = sexForColumn[sampleIdx];
            if ( sampleSex == 1 ) {
                maleCounts[rdCNVal] += 1;
            } else if ( sampleSex == 2 ) {
                femaleCounts[rdCNVal] += 1;
            }
        }
        final double maleMedian = calcMedian(maleCounts);
        double femaleMedian = calcMedian(femaleCounts);
        return maleMedian == 1. && (isY ? femaleMedian == 0. : femaleMedian == 2.);
    }

    // visible for testing
    static double calcMedian( final int[] counts ) {
        final double target = (counts[0] + counts[1] + counts[2] + counts[3]) / 2.;
        if ( target == 0. ) {
            return Double.NaN;
        }
        int total = 0;
        for ( int iii = 0; iii < 4; ++iii ) {
            total += counts[iii];
            if ( total == target ) {
                return iii + .5;
            } else if ( total > target ) {
                return (double)iii;
            }
        }
        throw new IllegalStateException("we should never reach this statement");
    }

    private static Set<ByteSequence> readNoisyEventsFile( final String filename ) {
        final Set<ByteSequence> noisyEvents = new HashSet<>();
        try {
            final BufferedReader neRdr =
                    new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
            String line;
            while ( (line = neRdr.readLine()) != null ) {
                noisyEvents.add(new ByteSequence(line));
            }
        } catch ( final IOException ioe ) {
            throw new RuntimeException("can't read noisy events file " + filename);
        }
        return noisyEvents;
    }

    private static Set<ByteSequence> readBothSidesFile( final String filename ) {
        final Set<ByteSequence> bothsidesEvents = new HashSet<>();
        try {
            final BufferedReader bsRdr =
                    new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
            String line;
            while ( (line = bsRdr.readLine()) != null ) {
                final String lastCol = line.substring(line.lastIndexOf('\t') + 1);
                bothsidesEvents.add(new ByteSequence(lastCol));
            }
        } catch ( final IOException ioe ) {
            throw new RuntimeException("can't read bothsides support file " + filename);
        }
        return bothsidesEvents;
    }

    private static int[] readPedFile( final String pedFilename, List<ByteSequence> sampleNames ) {
        final int nCols = sampleNames.size() - 9;
        final Map<ByteSequence, Integer> sexForSampleMap = new HashMap<>(2*nCols);
        final int[] sexForSample = new int[nCols];
        try {
            final BufferedReader pedRdr =
                    new BufferedReader(new InputStreamReader(new FileInputStream(pedFilename)));
            final Pattern tabPattern = Pattern.compile("\\t");
            String line;
            while ( (line = pedRdr.readLine()) != null ) {
                final Scanner scanner = new Scanner(line).useDelimiter(tabPattern);
                scanner.next(); // family ignored
                final String sampleName = scanner.next();
                scanner.next(); // mom ignored
                scanner.next(); // pop ignored
                final int sex = scanner.nextInt();
                sexForSampleMap.put(new ByteSequence(sampleName), sex);
            }
        } catch ( final IOException ioe ) {
            throw new RuntimeException("can't read " + pedFilename, ioe);
        }
        for ( int col = 0; col < nCols; ++col ) {
            final Integer sex = sexForSampleMap.get(sampleNames.get(col + 9));
            if ( sex == null ) {
                throw new RuntimeException("can't determine sex for sample " + sampleNames.get(col));
            }
            sexForSample[col] = sex;
        }
        return sexForSample;
    }
}
