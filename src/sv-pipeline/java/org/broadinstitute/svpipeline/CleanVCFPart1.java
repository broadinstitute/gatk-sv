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
    private static final ByteSequence UNRESOLVED_SUFFIX = new ByteSequence(";UNRESOLVED");
    private static final ByteSequence HIGH_SR_BACKGROUND_SUFFIX = new ByteSequence(";HIGH_SR_BACKGROUND");
    private static final ByteSequence DEL_VALUE = new ByteSequence("DEL");
    private static final ByteSequence DUP_VALUE = new ByteSequence("DUP");
    private static final ByteSequence RDCN_VALUE = new ByteSequence("RD_CN");
    private static final ByteSequence MISSING_VALUE = new ByteSequence(".");
    private static final ByteSequence MISSING_GENOTYPE = new ByteSequence("./.");
    private static final ByteSequence MISSING_SUFFIX = new ByteSequence(":.");
    private static final ByteSequence GT_REF_REF = new ByteSequence("0/0");
    private static final ByteSequence GT_REF_ALT = new ByteSequence("0/1");
    private static final ByteSequence GT_ALT_ALT = new ByteSequence("1/1");

    private static final int MIN_ALLOSOME_EVENT_SIZE = 5000;

    public static void main( final String[] args ) {
        if ( args.length != 5 ) {
            System.err.println("Usage: java CleanVCFPart1 INVCFFILE PEDIGREES XCHR YCHR NOISYEVENTS");
            System.exit(1);
        }
        final VCFParser parser = new VCFParser(args[0]);
        final ByteSequence xChrName = new ByteSequence(args[2]);
        final ByteSequence yChrName = new ByteSequence(args[3]);
        final Set<ByteSequence> noisyEvents = readNoisyEventsFile(args[4]);
        try ( final OutputStream os
                      = new BufferedOutputStream(new FileOutputStream(FileDescriptor.out)) ) {
            int[] sexForSample = null;
            while ( parser.hasMetadata() ) {
                final Metadata metadata = parser.nextMetaData();
                if ( metadata instanceof Columns ) {
                    final Columns cols = ((Columns)metadata);
                    sexForSample = readPedFile(args[1], cols.getValue());
                    os.write(("##FILTER=<ID=HIGH_SR_BACKGROUND,Description=\"High number of "
                            + "SR splits in background samples indicating messy region\">\n")
                                .getBytes(StandardCharsets.UTF_8));
                    os.write("##FILTER=<ID=UNRESOLVED,Description=\"Variant is unresolved\">\n"
                                .getBytes(StandardCharsets.UTF_8));
                } else if ( metadata instanceof KeyAttributes ) {
                    final KeyAttributes keyAttrs = (KeyAttributes)metadata;
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
                final int fmtIdx = record.getFormatIndex(EV_VALUE);
                if ( fmtIdx >= 0 ) {
                    final List<ByteSequence> genotypes = record.getGenotypes();
                    final int nGenotypes = genotypes.size();
                    for ( int gIdx = 0; gIdx < nGenotypes; ++gIdx ) {
                        final ByteSequence genotype = genotypes.get(gIdx);
                        final List<ByteSequence> genotypeVals = genotype.split(':');
                        final ByteSequence evVal = genotypeVals.get(fmtIdx);
                        genotypes.set(gIdx, genotype.replace(evVal, EV_VALS[evVal.asInt()]));
                    }
                }

                // move the SVTYPE to the ALT field (except for MEs)
                Map<ByteSequence, ByteSequence> infoMap = record.getInfoAsMap();
                final ByteSequence svType = infoMap.get(SVTYPE_KEY);
                if ( !record.getAlt().contains(ME_VALUE) ) {
                    if ( svType != null ) {
                        record.setAlt(new ByteSequence(LT_VALUE, svType, GT_VALUE));
                    }
                }
                record.setRef(N_VALUE);

                // move varGQ info field to quality column
                final ByteSequence varGQ = infoMap.get(VARGQ_KEY);
                if ( varGQ != null ) {
                    record.setQuality(varGQ);
                    record.removeInfoField(VARGQ_KEY);
                }

                // remove MULTIALLELIC and UNRESOLVED info flags
                if ( infoMap.containsKey(MULTIALLELIC_KEY) ) {
                    record.removeInfoField(MULTIALLELIC_KEY);
                }
                if ( infoMap.containsKey(UNRESOLVED_KEY) ) {
                    record.setFilter(new ByteSequence(record.getFilter(), UNRESOLVED_SUFFIX));
                    record.removeInfoField(UNRESOLVED_KEY);
                }

                // mark noisy events
                if ( noisyEvents.contains(record.getID()) ) {
                    record.setFilter(new ByteSequence(record.getFilter(), HIGH_SR_BACKGROUND_SUFFIX));
                }

                // fix genotypes on allosomes
                final boolean isY;
                if ( (isY = yChrName.equals(record.getChromosome())) ||
                        xChrName.equals(record.getChromosome())) {
                    final List<ByteSequence> genotypes = record.getGenotypes();
                    final int nSamples = genotypes.size();
                    final List<List<ByteSequence>> splitGenotypes = new ArrayList<>(nSamples);
                    for ( final ByteSequence genotype : genotypes ) {
                        splitGenotypes.add(genotype.split(':'));
                    }
                    final boolean isDel = DEL_VALUE.equals(svType);
                    final ByteSequence end = infoMap.get(END_KEY);
                    final int rdCNIndex = record.getFormatIndex(RDCN_VALUE);
                    boolean adjustMale = false;
                    if ( (isDel || DUP_VALUE.equals(svType)) && rdCNIndex >= 0 && end != null &&
                            end.asInt() + 1 - record.getPosition() > MIN_ALLOSOME_EVENT_SIZE ) {
                        adjustMale = isRevisableEvent(splitGenotypes, rdCNIndex, sexForSample, isY);
                    }
                    ByteSequence emptyGenotype = null;
                    for ( int sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx ) {
                        final int sampleSex = sexForSample[sampleIdx];
                        if ( sampleSex == 1 ) {
                            if ( adjustMale ) {
                                final List<ByteSequence> genotypeFields = splitGenotypes.get(sampleIdx);
                                final ByteSequence rdCN = genotypeFields.get(rdCNIndex);
                                if ( rdCN.equals(MISSING_VALUE) ) {
                                    continue;
                                }
                                final int rdCNVal = rdCN.asInt();
                                genotypeFields.set(rdCNIndex, new ByteSequence(Integer.toString(rdCNVal + 1)));
                                if ( isDel ) {
                                    if ( rdCNVal >= 1 ) genotypeFields.set(0, GT_REF_REF);
                                    else if ( rdCNVal == 0 ) genotypeFields.set(0, GT_REF_ALT);
                                } else {
                                    if ( rdCNVal <= 1 ) genotypeFields.set(0, GT_REF_REF);
                                    else if ( rdCNVal == 2 ) genotypeFields.set(0, GT_REF_ALT);
                                    else genotypeFields.set(0, GT_ALT_ALT);
                                }
                                genotypes.set(sampleIdx, new ByteSequence(genotypeFields, ':'));
                            }
                        } else if ( sampleSex == 2 ) {
                            if ( isY ) {
                                if ( emptyGenotype == null ) {
                                    emptyGenotype = MISSING_GENOTYPE;
                                    int nFields = splitGenotypes.get(sampleIdx).size();
                                    while ( --nFields > 0 ) {
                                        emptyGenotype = new ByteSequence(emptyGenotype, MISSING_SUFFIX);
                                    }
                                }
                                genotypes.set(sampleIdx, emptyGenotype);
                            }
                        } else {
                            final ByteSequence curGenotype = genotypes.get(sampleIdx);
                            final ByteSequence curGenotypeField =
                                    splitGenotypes.get(sampleIdx).get(0); // GT field is always first
                            final ByteSequence newGenotype =
                                    curGenotype.replace(curGenotypeField, MISSING_GENOTYPE);
                            genotypes.set(sampleIdx, newGenotype);
                        }
                    }
                }

                record.write(os);
            }
        } catch ( final IOException ioe ) {
            throw new RuntimeException("Can't write to stdout", ioe);
        }
    }

    private static boolean isRevisableEvent( final List<List<ByteSequence>> splitGenotypes,
                                             final int rdCNIndex,
                                             final int[] sexForColumn,
                                             final boolean isY ) {
        final int[] maleCounts = new int[4];
        final int[] femaleCounts = new int[4];
        final int nSamples = splitGenotypes.size();
        for ( int sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx ) {
            final ByteSequence rdCN = splitGenotypes.get(sampleIdx).get(rdCNIndex);
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
        final double maleTarget =
                (maleCounts[0] + maleCounts[1] + maleCounts[2] + maleCounts[3])/2.;
        double maleMedian = 0.;
        int counts = 0;
        for ( int iii = 0; iii < 4; ++iii ) {
            counts += maleCounts[iii];
            if ( counts == maleTarget ) {
                maleMedian = iii + .5;
                break;
            } else if ( counts > maleTarget ) {
                maleMedian = iii;
                break;
            }
        }
        final double femaleTarget =
                (femaleCounts[0] + femaleCounts[1] + femaleCounts[2] + femaleCounts[3])/2.;
        double femaleMedian = 0.;
        counts = 0;
        for ( int iii = 0; iii < 4; ++iii ) {
            counts += femaleCounts[iii];
            if ( counts == femaleTarget ) {
                femaleMedian = iii + .5;
                break;
            } else if ( counts > femaleTarget ) {
                femaleMedian = iii;
                break;
            }
        }
        return (!isY && maleMedian==1. && femaleMedian==2.) ||
                (isY && maleMedian==1. && femaleMedian==0.);
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
