/*
 * This file is part of DIA-Umpire.
 *
 * DIA-Umpire is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * DIA-Umpire is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with DIA-Umpire. If not, see <https://www.gnu.org/licenses/>.
 *
 */
package MSUmpire.SpectrumParser;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.*;
import java.util.stream.Stream;
import java.util.zip.Deflater;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.eclipse.collections.impl.list.mutable.primitive.FloatArrayList;
import org.eclipse.collections.impl.list.mutable.primitive.IntArrayList;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scan.props.PrecursorInfo;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.filetypes.AbstractLCMSDataSource;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;
import umich.ms.fileio.filetypes.thermo.ThermoRawFile;
import umich.ms.msfiletoolbox.MsftbxInfo;

/* * 
 * mzXML parser
 */
/**
 *
 * @author Chih-Chiang Tsou
 */
public final class mzXMLParser  extends SpectrumParserBase{
    public static void main(final String[] args) throws Exception{
//        to_mzXML("/home/ci/tmp/tkQE170512_U_ThermoFixed_DIA_01_Q3.mzXML");
        {
            final long a = System.currentTimeMillis();
            to_mzXML("/home/ci/DIA-U_batmass_io_test/JHU_LM_DIA_Pancreatic_2A6_02.mzML", 10);
            System.out.println("time taken input to mzML:" + (System.currentTimeMillis() - a) / 1000. + "s");
        }
        try (final BufferedWriter bw = Files.newBufferedWriter(Paths.get("/home/ci/DIA-U_batmass_io_test/output/JHU_LM_DIA_Pancreatic_2A6_02_Q1_test.mzML"), StandardCharsets.US_ASCII)) {
            final long a = System.currentTimeMillis();
            to_mzML(Paths.get("/home/ci/DIA-U_batmass_io_test/output/JHU_LM_DIA_Pancreatic_2A6_02_Q1.mgf"), bw, 10);
            System.out.println("time taken mgf to mzML:"+(System.currentTimeMillis()-a)/1000.+"s");
            test_to_mzML("/home/ci/DIA-U_batmass_io_test/output/JHU_LM_DIA_Pancreatic_2A6_02_Q1_test.mzML");
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }

    }

    public static void test_to_mzML(final String path) {
        try {
            final int numThreads = 10;
            final int parsingTimeout = 4;
            ScanCollectionDefault scans = new ScanCollectionDefault();
            scans.setDefaultStorageStrategy(StorageStrategy.STRONG);
            scans.isAutoloadSpectra(true);

            final MZMLFile source = new MZMLFile(path);

            source.setExcludeEmptyScans(true);
            source.setNumThreadsForParsing(numThreads);
            source.setParsingTimeout(parsingTimeout);
            scans.setDataSource(source);

            scans.loadData(LCMSDataSubset.WHOLE_RUN);
            final TreeMap<Integer, IScan> num2scan = scans.getMapNum2scan();
            System.out.println("num2scan.size() = " + num2scan.size());
//            System.out.println("num2scan = " + num2scan);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void to_mzML(final Path mgf, final BufferedWriter xmloutput, final int numThreads) throws IOException, java.security.NoSuchAlgorithmException {
        final String head_format_str = "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" +
                "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">\n" +
                "  <mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" id=\"%s\" version=\"1.1.0\">\n" +
                "    <cvList count=\"2\">\n" +
                "      <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"3.60.0\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\"/>\n" +
                "      <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"12:10:2011\" URI=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\"/>\n" +
                "    </cvList>\n" +
                "    <fileDescription>\n" +
                "      <fileContent>\n" +
                "        <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n" +
                "        <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n" +
                "      </fileContent>\n" +
                "    </fileDescription>\n" +
                "    <softwareList count=\"1\">\n" +
                "      <software id=\"pwiz_3.0.6002_x0020__x0028_TPP_x0020_v5.0.0_x0020_Typhoon_x002c__x0020_Build_x0020_201610181405-exported_x0020__x0028_Linux-x86_64_x0029__x0029_\" version=\"3.0.6002 (TPP v5.0.0 Typhoon, Build 201610181405-exported (Linux-x86_64))\">\n" +
                "        <cvParam cvRef=\"MS\" accession=\"MS:1000615\" name=\"ProteoWizard software\" value=\"\"/>\n" +
                "      </software>\n" +
                "    </softwareList>\n" +
                "    <instrumentConfigurationList count=\"1\">\n" +
                "      <instrumentConfiguration id=\"IC\">\n" +
                "        <cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\" value=\"\"/>\n" +
                "      </instrumentConfiguration>\n" +
                "    </instrumentConfigurationList>\n" +
                "    <dataProcessingList count=\"1\">\n" +
                "      <dataProcessing id=\"pwiz_Reader_conversion\">\n" +
                "        <processingMethod order=\"0\" softwareRef=\"pwiz_3.0.6002_x0020__x0028_TPP_x0020_v5.0.0_x0020_Typhoon_x002c__x0020_Build_x0020_201610181405-exported_x0020__x0028_Linux-x86_64_x0029__x0029_\">\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\" value=\"\"/>\n" +
                "        </processingMethod>\n" +
                "      </dataProcessing>\n" +
                "    </dataProcessingList>\n" +
                "    <run id=\"%s\" defaultInstrumentConfigurationRef=\"IC\">\n" +
                "      <spectrumList count=\"%d\" defaultDataProcessingRef=\"pwiz_Reader_conversion\">\n";
        final String spectrum_indent = "        ";
        final String spectrum_format_str = "<spectrum index=\"%d\" id=\"controllerType=0 controllerNumber=1 scan=%d\" defaultArrayLength=\"%d\">\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000796\" name=\"spectrum title\" value=\"%s\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\" value=\"\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000528\" name=\"lowest observed m/z\" value=\"%f\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000527\" name=\"highest observed m/z\" value=\"%f\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000285\" name=\"total ion current\" value=\"%f\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000504\" name=\"base peak m/z\" value=\"%f\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000505\" name=\"base peak intensity\" value=\"%f\"/>\n" +
                "          <scanList count=\"1\">\n" +
                "            <cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>\n" +
                "            <scan>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"%f\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>\n" +
                "            </scan>\n" +
                "          </scanList>\n" +
                "          <precursorList count=\"1\">\n" +
                "            <precursor>\n" +
                "              <selectedIonList count=\"1\">\n" +
                "                <selectedIon>\n" +
                "                  <cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"%f\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n" +
                "                  <cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"%d\"/>\n" +
                "                </selectedIon>\n" +
                "              </selectedIonList>\n" +
                "              <activation>\n" +
                "              </activation>\n" +
                "            </precursor>\n" +
                "          </precursorList>\n" +
                "          <binaryDataArrayList count=\"2\">\n" +
                "            <binaryDataArray encodedLength=\"%d\">\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n" +
                "              <binary>%s</binary>\n" +
                "            </binaryDataArray>\n" +
                "            <binaryDataArray encodedLength=\"%d\">\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n" +
                "              <binary>%s</binary>\n" +
                "            </binaryDataArray>\n" +
                "          </binaryDataArrayList>\n" +
                "        </spectrum>\n";

        final MessageDigest sha1 = MessageDigest.getInstance("SHA-1");
        final String run_id0 = mgf.getFileName().toString();
        final String run_id = run_id0.substring(0, run_id0.length() - 4);
        final List<Future<String>> futures = new ArrayList<>();
        final ForkJoinPool fjp = new ForkJoinPool(numThreads);

        try (final BufferedReader br = Files.newBufferedReader(mgf)) {
            String line;
            List<String> ions_entry = null;
            while ((line = br.readLine()) != null) {
                if (line.equals("BEGIN IONS")) {
                    ions_entry = new ArrayList<>();
                    line = br.readLine();
                    if (line == null)
                        break;
                }
                if (!line.equals("END IONS")) {
                    if (ions_entry != null)
                        ions_entry.add(line);
                } else {
                    final List<String> l = ions_entry;
                    final int index = futures.size();
                    futures.add(fjp.submit(()-> {
                        final ListIterator<String> li = l.listIterator();
                        final double pepmass = Double.parseDouble(li.next().split("=")[1]);
                        final int charge = Integer.parseInt(li.next().split("[=+]")[1]);
                        final double rtinseconds = Double.parseDouble(li.next().split("=")[1]);
                        final String title = li.next().split("=")[1];
                        final int defaultArrayLength = l.size() - li.nextIndex();
                        final FloatArrayList mzarr = new FloatArrayList(defaultArrayLength);
                        final FloatArrayList intensityarr = new FloatArrayList(defaultArrayLength);
                        while(li.hasNext()){
                            final String[] mz_int = li.next().split(" ");
                            mzarr.add(Float.parseFloat(mz_int[0]));
                            intensityarr.add(Float.parseFloat(mz_int[1]));
                        }

                        final ByteBuffer bb = ByteBuffer.allocate(defaultArrayLength * Float.BYTES);
                        bb.order(ByteOrder.LITTLE_ENDIAN);
                        mzarr.forEach(bb::putFloat);
                        final String base64_mz_array = Base64.getEncoder().encodeToString(bb.array());
                        bb.clear();
                        intensityarr.forEach(bb::putFloat);
                        final String base64_intensity_array = Base64.getEncoder().encodeToString(bb.array());
                        final float base_peak_intensity = intensityarr.max();
                        final String spectrum_xml = spectrum_indent + String.format(spectrum_format_str,
                                index, index + 1, defaultArrayLength,
                                title,
                                mzarr.min(), mzarr.max(),
                                intensityarr.sum(),
                                mzarr.get(intensityarr.indexOf(base_peak_intensity)), base_peak_intensity,
                                rtinseconds, pepmass, charge,
                                base64_mz_array.length(), base64_mz_array,
                                base64_intensity_array.length(), base64_intensity_array);
                        return spectrum_xml;
                    }));
                    ions_entry = null;
                }
            }
        }

        int char_count = 0;
        final int spectrumList_count = futures.size();
        final long[] index_spectrum_offset = new long[spectrumList_count];
        final String head_xml = String.format(head_format_str, run_id, run_id, spectrumList_count);
        sha1.update(head_xml.getBytes(StandardCharsets.US_ASCII));
        char_count += head_xml.length();
        xmloutput.write(head_xml);

        if(!true)
        for (int index = 0; index < futures.size(); ++index) {
            final String spectrum_xml;
            try {
                spectrum_xml = futures.get(index).get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException(e);
            }
            index_spectrum_offset[index] = char_count + spectrum_indent.length();
            sha1.update(spectrum_xml.getBytes(StandardCharsets.US_ASCII));
            char_count += spectrum_xml.length();
            xmloutput.write(spectrum_xml);
        }
        final ForkJoinTask<?> fjt_write = fjp.submit(() -> {
            for (final Future<String> future : futures) {
                final String spectrum_xml;
                try {
                    spectrum_xml = future.get();
                } catch (InterruptedException | ExecutionException e) {
                    throw new RuntimeException(e);
                }
                try {
                    xmloutput.write(spectrum_xml);
                } catch (IOException e) {
                    throw new UncheckedIOException(e);
                }
            }
        });
        for (int index = 0; index < futures.size(); ++index) {
            final String spectrum_xml;
            try {
                spectrum_xml = futures.get(index).get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException(e);
            }
            index_spectrum_offset[index] = char_count + spectrum_indent.length();
            sha1.update(spectrum_xml.getBytes(StandardCharsets.US_ASCII));
            char_count += spectrum_xml.length();
        }

        try {
            fjt_write.get();
        } catch (InterruptedException | ExecutionException e) {
            throw new RuntimeException(e);
        }

        {
            final String tmp1 = "      </spectrumList>\n" +
                    "    </run>\n" +
                    "  </mzML>\n" +
                    "  ";
            sha1.update(tmp1.getBytes(StandardCharsets.US_ASCII));
            char_count += tmp1.length();
            xmloutput.write(tmp1);
        }
        final long indexListOffset = char_count;
        {
            final String tmp2 = "<indexList count=\"2\">\n" +
                    "    <index name=\"spectrum\">\n";
            sha1.update(tmp2.getBytes(StandardCharsets.US_ASCII));
            xmloutput.write(tmp2);
        }

        final StringBuilder sb = new StringBuilder();
        for (int index = 0; index < spectrumList_count; index++) {
            final String offset_xml = String.format("      <offset idRef=\"controllerType=0 controllerNumber=1 scan=%d\">%d</offset>\n",
                    index + 1, index_spectrum_offset[index]);
            sb.append(offset_xml);
        }
        sb.append(String.format("    </index>\n" +
                "    <index name=\"chromatogram\">\n" +
                "    </index>\n" +
                "  </indexList>\n" +
                "  <indexListOffset>%d</indexListOffset>\n" +
                "  <fileChecksum>", indexListOffset));

        sha1.update(sb.toString().getBytes(StandardCharsets.US_ASCII));
        sb.append(String.format("%s</fileChecksum>\n" +
                        "</indexedmzML>\n",
//                new String(org.apache.commons.codec.binary.Hex.encodeHex(sha1.digest()))));
                javax.xml.bind.DatatypeConverter.printHexBinary(sha1.digest()).toLowerCase(Locale.ROOT)));
        xmloutput.write(sb.toString());
        xmloutput.close();
    }

    public static void to_mzML_sequential(final Path mgf, BufferedWriter xmloutput) throws IOException, java.security.NoSuchAlgorithmException {
        final String head_format_str = "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" +
                "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">\n" +
                "  <mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" id=\"%s\" version=\"1.1.0\">\n" +
                "    <cvList count=\"2\">\n" +
                "      <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"3.60.0\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\"/>\n" +
                "      <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"12:10:2011\" URI=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\"/>\n" +
                "    </cvList>\n" +
                "    <fileDescription>\n" +
                "      <fileContent>\n" +
                "        <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n" +
                "        <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n" +
                "      </fileContent>\n" +
                "    </fileDescription>\n" +
                "    <softwareList count=\"1\">\n" +
                "      <software id=\"pwiz_3.0.6002_x0020__x0028_TPP_x0020_v5.0.0_x0020_Typhoon_x002c__x0020_Build_x0020_201610181405-exported_x0020__x0028_Linux-x86_64_x0029__x0029_\" version=\"3.0.6002 (TPP v5.0.0 Typhoon, Build 201610181405-exported (Linux-x86_64))\">\n" +
                "        <cvParam cvRef=\"MS\" accession=\"MS:1000615\" name=\"ProteoWizard software\" value=\"\"/>\n" +
                "      </software>\n" +
                "    </softwareList>\n" +
                "    <instrumentConfigurationList count=\"1\">\n" +
                "      <instrumentConfiguration id=\"IC\">\n" +
                "        <cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\" value=\"\"/>\n" +
                "      </instrumentConfiguration>\n" +
                "    </instrumentConfigurationList>\n" +
                "    <dataProcessingList count=\"1\">\n" +
                "      <dataProcessing id=\"pwiz_Reader_conversion\">\n" +
                "        <processingMethod order=\"0\" softwareRef=\"pwiz_3.0.6002_x0020__x0028_TPP_x0020_v5.0.0_x0020_Typhoon_x002c__x0020_Build_x0020_201610181405-exported_x0020__x0028_Linux-x86_64_x0029__x0029_\">\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\" value=\"\"/>\n" +
                "        </processingMethod>\n" +
                "      </dataProcessing>\n" +
                "    </dataProcessingList>\n" +
                "    <run id=\"%s\" defaultInstrumentConfigurationRef=\"IC\">\n" +
                "      <spectrumList count=\"%d\" defaultDataProcessingRef=\"pwiz_Reader_conversion\">\n";
        final String spectrum_indent = "        ";
        final String spectrum_format_str = "<spectrum index=\"%d\" id=\"controllerType=0 controllerNumber=1 scan=%d\" defaultArrayLength=\"%d\">\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000796\" name=\"spectrum title\" value=\"%s\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\" value=\"\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000528\" name=\"lowest observed m/z\" value=\"%f\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000527\" name=\"highest observed m/z\" value=\"%f\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000285\" name=\"total ion current\" value=\"%f\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000504\" name=\"base peak m/z\" value=\"%f\"/>\n" +
                "          <cvParam cvRef=\"MS\" accession=\"MS:1000505\" name=\"base peak intensity\" value=\"%f\"/>\n" +
                "          <scanList count=\"1\">\n" +
                "            <cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>\n" +
                "            <scan>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"%f\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>\n" +
                "            </scan>\n" +
                "          </scanList>\n" +
                "          <precursorList count=\"1\">\n" +
                "            <precursor>\n" +
                "              <selectedIonList count=\"1\">\n" +
                "                <selectedIon>\n" +
                "                  <cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"%f\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n" +
                "                  <cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"%d\"/>\n" +
                "                </selectedIon>\n" +
                "              </selectedIonList>\n" +
                "              <activation>\n" +
                "              </activation>\n" +
                "            </precursor>\n" +
                "          </precursorList>\n" +
                "          <binaryDataArrayList count=\"2\">\n" +
                "            <binaryDataArray encodedLength=\"%d\">\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n" +
                "              <binary>%s</binary>\n" +
                "            </binaryDataArray>\n" +
                "            <binaryDataArray encodedLength=\"%d\">\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>\n" +
                "              <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n" +
                "              <binary>%s</binary>\n" +
                "            </binaryDataArray>\n" +
                "          </binaryDataArrayList>\n" +
                "        </spectrum>\n";

        final MessageDigest sha1 = MessageDigest.getInstance("SHA-1");
        final String run_id0 = mgf.getFileName().toString();
        final String run_id = run_id0.substring(0, run_id0.length() - 4);
        final int spectrumList_count;
        try (final BufferedReader br = Files.newBufferedReader(mgf)) {
            int count = 0;
            String line;
            while ((line = br.readLine()) != null)
                if (line.equals("BEGIN IONS"))
                    ++count;
            spectrumList_count = count;
        }
        long char_count = 0;
        final long[] index_spectrum_offset = new long[spectrumList_count];

        final String head_xml = String.format(head_format_str, run_id, run_id, spectrumList_count);
        sha1.update(head_xml.getBytes(StandardCharsets.US_ASCII));
        char_count += head_xml.length();
        xmloutput.write(head_xml);

        final BufferedReader br = Files.newBufferedReader(mgf);
        String line;

        for (int index = 0; index < spectrumList_count; index++) {

            double pepmass = Double.NaN;
            int charge;
            double rtinseconds;
            String title;

            line = br.readLine();
            if (line == null)
                break;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("PEPMASS")) {
                    pepmass = Double.parseDouble(line.split("=")[1]);
                    break;
                }
            }
            charge = Integer.parseInt(br.readLine().split("[=+]")[1]);
            rtinseconds = Double.parseDouble(br.readLine().split("=")[1]);
            title = br.readLine().split("=")[1];
            final FloatArrayList mzarr = new FloatArrayList();
            final FloatArrayList intensityarr = new FloatArrayList();
            while (!(line = br.readLine()).equals("END IONS")) {
                final String[] mz_int = line.split(" ");
                mzarr.add(Float.parseFloat(mz_int[0]));
                intensityarr.add(Float.parseFloat(mz_int[1]));
            }
            final int defaultArrayLength = mzarr.size();

            final ByteBuffer bb = ByteBuffer.allocate(defaultArrayLength * Float.BYTES);
            bb.order(ByteOrder.LITTLE_ENDIAN);
            mzarr.forEach(bb::putFloat);
            final String base64_mz_array = Base64.getEncoder().encodeToString(bb.array());
            bb.clear();
            intensityarr.forEach(bb::putFloat);
            final String base64_intensity_array = Base64.getEncoder().encodeToString(bb.array());

            final String spectrum_xml = spectrum_indent + String.format(spectrum_format_str,
                    index, index + 1, defaultArrayLength,
                    title,
                    mzarr.min(), mzarr.max(),
                    intensityarr.sum(),
                    mzarr.get(intensityarr.indexOf(intensityarr.max())), intensityarr.max(),
                    rtinseconds, pepmass, charge,
                    base64_mz_array.length(), base64_mz_array,
                    base64_intensity_array.length(), base64_intensity_array);

            index_spectrum_offset[index] = char_count + spectrum_indent.length();
            sha1.update(spectrum_xml.getBytes(StandardCharsets.US_ASCII));
            char_count += spectrum_xml.length();
            xmloutput.write(spectrum_xml);
        }
        {
            final String tmp1 = "      </spectrumList>\n" +
                    "    </run>\n" +
                    "  </mzML>\n" +
                    "  ";
            sha1.update(tmp1.getBytes(StandardCharsets.US_ASCII));
            char_count += tmp1.length();
            xmloutput.write(tmp1);
        }
        final long indexListOffset = char_count;
        {
            final String tmp2 = "<indexList count=\"2\">\n" +
                    "    <index name=\"spectrum\">\n";
            sha1.update(tmp2.getBytes(StandardCharsets.US_ASCII));
            xmloutput.write(tmp2);
        }

        final StringBuilder sb = new StringBuilder();
        for (int index = 0; index < spectrumList_count; index++) {
            final String offset_xml = String.format("      <offset idRef=\"controllerType=0 controllerNumber=1 scan=%d\">%d</offset>\n",
                    index + 1, index_spectrum_offset[index]);
            sb.append(offset_xml);
        }
        sb.append(String.format("    </index>\n" +
                "    <index name=\"chromatogram\">\n" +
                "    </index>\n" +
                "  </indexList>\n" +
                "  <indexListOffset>%d</indexListOffset>\n" +
                "  <fileChecksum>", indexListOffset));

        sha1.update(sb.toString().getBytes(StandardCharsets.US_ASCII));
        sb.append(String.format("%s</fileChecksum>\n" +
                        "</indexedmzML>\n",
//                new String(org.apache.commons.codec.binary.Hex.encodeHex(sha1.digest()))));
                javax.xml.bind.DatatypeConverter.printHexBinary(sha1.digest()).toLowerCase(Locale.ROOT)));
        xmloutput.write(sb.toString());
        xmloutput.close();
    }

    public static String to_mzXML(final String path, final int numThreads) throws Exception {
        final int parsingTimeout = 4;
        ScanCollectionDefault scans = new ScanCollectionDefault();
        scans.setDefaultStorageStrategy(StorageStrategy.STRONG);
        scans.isAutoloadSpectra(true);

        final String basename = path.substring(0, path.lastIndexOf("."));
        final String ext = path.substring(path.lastIndexOf(".") + 1);
        final AbstractLCMSDataSource<?> source;
        switch (ext.toLowerCase()) {
            case "mzml":
                source = new MZMLFile(path);
                break;
            case "mzxml":
                source = new MZXMLFile(path);
                break;
            case "raw":
                System.out.println("Batmass-IO version " + MsftbxInfo.getVersion());
                source = new ThermoRawFile(path);
                break;
            default:
                throw new RuntimeException("Unrecognized file extension: " + ext);
        }

        source.setExcludeEmptyScans(true);
        source.setNumThreadsForParsing(numThreads);
        source.setParsingTimeout(parsingTimeout);
        scans.setDataSource(source);

        scans.loadData(LCMSDataSubset.WHOLE_RUN);
        final TreeMap<Integer, IScan> num2scan = scans.getMapNum2scan();
        final Path ret = Paths.get(basename + ".mzXML");
        Files.deleteIfExists(ret);
        ret.toFile().deleteOnExit();
        final OutputStreamWriter writer = new OutputStreamWriter(Files.newOutputStream(ret, StandardOpenOption.CREATE_NEW), StandardCharsets.ISO_8859_1);
        writer.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" +
                "<mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.2\"\n" +
                "       xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n" +
                "       xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.2 http://sashimi.sourceforge.net/schema_revision/mzXML_3.2/mzXML_idx_3.2.xsd\">\n" +
                "  <msRun>\n");
//        System.out.println("num2scan.size() = " + num2scan.size());
        final Stream<String> stringStream = num2scan.entrySet().parallelStream().map(entry -> {
            final StringBuilder sb = new StringBuilder();
            final IScan value = entry.getValue();
            final ISpectrum spectrum = value.getSpectrum();
            sb.append(String.format("    <scan num=\"%d\"\n" +
                            "          scanType=\"%s\"\n" +
                            "          centroided=\"%d\"\n" +
                            "          msLevel=\"%d\"\n" +
                            "          peaksCount=\"%d\"\n" +
                            "          polarity=\"%s\"\n" +
                            "          retentionTime=\"PT%s\"\n" +
                            "          lowMz=\"%f\"\n" +
                            "          highMz=\"%f\"\n" +
                            "          basePeakMz=\"%f\"\n" +
                            "          basePeakIntensity=\"%f\"\n" +
                            "          totIonCurrent=\"%f\">\n",
                    value.getNum(),
                    value.getScanType() == null ? "Full" : value.getScanType().name(),
                    value.isCentroided() ? 1 : 0,
                    value.getMsLevel(),
                    spectrum.getMZs().length,
                    value.getPolarity().toString(),
                    value.getRt() * 60 + "S",
                    spectrum.getMinMZ(),
                    spectrum.getMaxMZ(),
                    value.getBasePeakMz(),
                    value.getBasePeakIntensity(),
                    value.getTic()));
            final PrecursorInfo precursor = value.getPrecursor();
            if (precursor != null) {
                final Double windowWideness = precursor.getMzRangeEnd() == null && precursor.getMzRangeStart() == null ?
                        null :
                        precursor.getMzRangeEnd() - precursor.getMzRangeStart();
                sb.append(String.format("      <precursorMz precursorIntensity=\"%f\"%s%s%s%s>%f</precursorMz>\n",
                        precursor.getIntensity() == null ? 0 : precursor.getIntensity(),
                        precursor.getParentScanNum() == null ? "" : String.format(" precursorScanNum=\"%d\"", precursor.getParentScanNum()),
                        precursor.getActivationInfo().getActivationMethod() == null ? "" : String.format(" activationMethod=\"%s\"", precursor.getActivationInfo().getActivationMethod().name),
                        precursor.getCharge() == null ? "" : String.format(" precursorCharge=\"%d\"", precursor.getCharge()),
                        windowWideness == null ? "" : String.format(" windowWideness=\"%f\"", windowWideness),
                        precursor.getMzTargetMono()));
            }
            final ByteBuffer byteBuffer = ByteBuffer.allocate(spectrum.getMZs().length * 2 * Float.BYTES);
            for (int i = 0; i < spectrum.getMZs().length; i++) {
                byteBuffer.putInt(Float.floatToRawIntBits((float) spectrum.getMZs()[i]));
                byteBuffer.putInt(Float.floatToRawIntBits((float) spectrum.getIntensities()[i]));
            }

            // Compress the bytes
            final byte[] output = new byte[byteBuffer.array().length << 1];
            final Deflater compresser = new Deflater();
            compresser.setInput(byteBuffer.array());
            compresser.finish();
            final int compressedDataLength = compresser.deflate(output);
            compresser.end();
            sb.append(String.format("      <peaks compressionType=\"zlib\"\n" +
                            "             compressedLen=\"%d\"\n" +
                            "             precision=\"32\"\n" +
                            "             byteOrder=\"network\"\n" +
                            "             contentType=\"m/z-int\">%s</peaks>\n",
                    compressedDataLength,
                    Base64.getEncoder().encodeToString(Arrays.copyOf(output, compressedDataLength))
            ));
            sb.append("    </scan>\n");
            return sb.toString();
        });
        stringStream.forEachOrdered(str -> {
            try {
                writer.write(str);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        });
        writer.write("  </msRun>\n");
        writer.write("</mzXML>");
        writer.close();
        return ret.toString();
    }
    public TreeMap<Integer, Long> ScanIndex=null;
    public mzXMLParser(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, DIA_Setting dIA_Setting, int NoCPUs) throws Exception {
        super(filename,parameter,datatype,dIA_Setting,NoCPUs);
        ReadElutionAndScanIndex();
    }

    //Parser elution time index and scan index and save them as binary files
    private void ReadElutionAndScanIndex() throws Exception {
        if (!FSScanPosRead()) {
            InferOffsets();
//            ParseIndex();
            WriteIndexSerialization();
        }
        if (!FSElutionIndexRead()) {
            ParseElutionIndex();
            FSElutionIndexWrite();
        }
    }

     //Wirte seralization file for scan index
    private void WriteIndexSerialization() {
        FSScanPosWrite();
    }
    
     //Wirte seralization file for scan index
    private void FSScanPosWrite() {
        try {
            LogManager.getRootLogger().debug("Writing ScanPos to file:" + FilenameUtils.removeExtension(filename) + ".ScanPosFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".ScanPosFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(ScanIndex);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    
    //Read seralization file for scan index
    private boolean FSScanPosRead() {
        if (!new File(FilenameUtils.removeExtension(filename) + ".ScanPosFS").exists()) {
            return false;
        }
        try {
            LogManager.getRootLogger().debug("Reading ScanPos:" + FilenameUtils.removeExtension(filename) + ".ScanPosFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".ScanPosFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            ScanIndex = (TreeMap<Integer, Long>) in.readObject();
            TotalScan = ScanIndex.size();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            LogManager.getRootLogger().debug("ScanIndex serialization file failed");
            //LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    /**
     * Infer scan tag offsets.
     * @throws IOException
     */
    private void InferOffsets() throws IOException {
        this.ScanIndex = new TreeMap<>();
        try (final InputStream is = Files.newInputStream(Paths.get(this.filename));) {
            final LongArrayList offsets = Infer_mzXML_offsets.infer_mzXML_scanno(is);
            this.TotalScan = offsets.size();
            for (int i=0;i<offsets.size()-1;++i){
                this.ScanIndex.put(i+1,offsets.get(i));
            }
            this.ScanIndex.put(Integer.MAX_VALUE, offsets.getLast());
        }
    }

    //Parse scan index at the bottom of mzXML file
    private void ParseIndex() throws FileNotFoundException, IOException {
        TotalScan = 0;
        ScanIndex = new TreeMap<>();
        try (RandomAccessFile fileHandler = new RandomAccessFile(filename, "r")) {
            StringBuilder sb = new StringBuilder();

            String CurrentLine = "";
            long currentLastPt = fileHandler.length() - 1;
            boolean indexexist=false;
            int linecount=0;
            while (!(CurrentLine.trim().startsWith("<index name=") | CurrentLine.trim().startsWith("</msRun>"))) {
                //Read backward
                for (long filePointer = currentLastPt; filePointer != -1; filePointer--) {
                    fileHandler.seek(filePointer);
                    int readByte = fileHandler.readByte();
                    if (readByte == 0xA) {
                        if (filePointer == currentLastPt) {
                            continue;
                        } else {
                            currentLastPt = filePointer;
                            break;
                        }
                    } else if (readByte == 0xD) {
                        if (filePointer == currentLastPt - 1) {
                            continue;
                        } else {
                            currentLastPt = filePointer;
                            break;
                        }
                    }
                    sb.append((char) readByte);
                }
                linecount++;
                CurrentLine = sb.reverse().toString();
                sb = new StringBuilder();

                if (CurrentLine.trim().startsWith("</index>")) {
                    indexexist = true;
                }

                if (!indexexist && linecount > 10) {
                    fileHandler.close();
                    LogManager.getRootLogger().debug("File : " + filename + " doesn't have index. the processing will stop.");
                    System.exit(1);
                }

                if (CurrentLine.trim().startsWith("<offset id")) {
                    int scanNo = Integer.parseInt(CurrentLine.substring(CurrentLine.indexOf("<offset id=\"") + 12).split("\"")[0]);
                    long index = (long) Long.parseLong(CurrentLine.substring(CurrentLine.indexOf(">") + 1, CurrentLine.indexOf("</offset>")));
                    if (index < 0) {
                        index = index + 2147483647l + 2147483648l;
                    }
                    if (ScanIndex.containsKey(scanNo + 1) && ScanIndex.get(scanNo + 1) == index) {
                        LogManager.getRootLogger().debug("File : " + filename + " index is not correct, ScanNo:" + scanNo + " and " + scanNo + 1 + " have same index");
                        LogManager.getRootLogger().debug("Please use indexmzXML from  TPP package to fix incorrect index of the mzXML file.");
                        LogManager.getRootLogger().debug("command: indexmzXML filename.mzXML");
                        System.exit(1);
                    }
                    ScanIndex.put(scanNo, index);
                } else if (CurrentLine.trim().startsWith("<indexOffset>")) {
                    long IndexEnd = (long) Long.parseLong(CurrentLine.substring(CurrentLine.indexOf("<indexOffset>") + 13, CurrentLine.indexOf("</indexOffset>")));
                    if (IndexEnd < 0) {
                        IndexEnd = IndexEnd + 2147483647l + 2147483648l;
                    }
                    ScanIndex.put(Integer.MAX_VALUE, IndexEnd);
                }
            }
            TotalScan = ScanIndex.size();
            sb = null;
            fileHandler.close();
        }
    }

     //Parse elution time-scan number mapping
    //For DIA data, isolation window ranges are parsed in this method
    private void ParseElutionIndex() throws Exception {
        
        if(ScanIndex==null | ScanIndex.isEmpty()){
            return;
        }

        try (RandomAccessFile fileHandler = new RandomAccessFile(filename, "r")) {
            Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();
            Long currentIdx = iter.next().getValue();
            while (iter.hasNext()) {
                long startposition = currentIdx;
                long nexposition = iter.next().getValue();
                currentIdx = nexposition;
                fileHandler.seek(startposition);
                
                byte[] bufr = new byte[(int) (nexposition - startposition)];
                fileHandler.readFully(bufr, 0, (int) (nexposition - startposition));

                String temp = new String(bufr);

                float rt = 0f;
                int scanno = 0;
                int mslevel = 0;
                //float precursorF=0f;
                if (!temp.contains("<scan")) {
                    fileHandler.close();
                    return;
                }

                if (temp.contains("<scan num=") && (temp.contains("retentionTime=\"PT"))) {
                    String substr = temp.substring(temp.indexOf("<scan num=") + 11);
                    scanno = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    
                    rt = Float.parseFloat(temp.substring(temp.indexOf("retentionTime=\"PT") + 17).split("S\"")[0]);
                    rt=rt/60f;
                    mslevel = Integer.parseInt(temp.substring(temp.indexOf("msLevel=") + 9, temp.indexOf("msLevel=") + 10));
                    if (temp.contains("scanType=\"calibration\"")) {
                        mslevel = -1;
                    }
                    if (mslevel == 1) {
                        NoMS1Scans++;                        
                        if (temp.contains("scanType=\"SIM\"") && datatype == SpectralDataType.DataType.WiSIM) {
                            int startidx = temp.indexOf("lowMz=\"") + 7;
                            int stopidx = startidx + 1;
                            for (int i = startidx + 1; i < temp.length(); i++) {
                                if (temp.charAt(i) == '\"') {
                                    stopidx = i;
                                    break;
                                }
                            }
                            float lowmz = Float.parseFloat(temp.substring(startidx, stopidx));
                            startidx = temp.indexOf("highMz=\"") + 8;
                            stopidx = startidx + 1;
                            for (int i = startidx + 1; i < temp.length(); i++) {
                                if (temp.charAt(i) == '\"') {
                                    stopidx = i;
                                    break;
                                }
                            }
                            float highmz = Float.parseFloat(temp.substring(startidx, stopidx));
                            for (XYData MS1win : dIA_Setting.MS1Windows.keySet()) {
                                if (MS1win.getX() <= lowmz && MS1win.getY() >= highmz) {
                                    dIA_Setting.MS1Windows.get(MS1win).add(scanno);
                                }
                            }
                        }                        
                    }
                    //If it is DIA data, parse isolation window ranges 
                    if (datatype != SpectralDataType.DataType.DDA) {
                        if (mslevel == 2) {
                            if (datatype == SpectralDataType.DataType.MSX) {
                                substr = temp;
                                while (substr.contains("</precursorMz>")) {
                                    int stopidx = substr.indexOf("</precursorMz>");
                                    int startidx = 0;
                                    for (int i = stopidx; i > 0; i--) {
                                        if (substr.charAt(i) == '>') {
                                            startidx = i + 1;
                                            break;
                                        }
                                    }
                                    float precursormz = Float.parseFloat(substr.substring(startidx, stopidx));

                                    startidx = substr.indexOf("windowWideness=\"") + 16;
                                    stopidx = startidx + 1;
                                    for (int i = startidx + 1; i < substr.length(); i++) {
                                        if (substr.charAt(i) == '\"') {
                                            stopidx = i;
                                            break;
                                        }
                                    }
                                    float windowwideness = Float.parseFloat(substr.substring(startidx, stopidx));                                    
                                    //Assuming the precursor m/z is at the center of isolation window, it's for Thermo MSX data
                                    float Loffset = windowwideness / 2f;
                                    float Roffset = windowwideness / 2f;

                                    if (!dIA_Setting.DIAWindows.containsKey(new XYData(precursormz - Loffset, precursormz + Roffset))) {
                                        ArrayList<Integer> scanList = new ArrayList<>();
                                        dIA_Setting.DIAWindows.put(new XYData(precursormz - Loffset, precursormz + Roffset), scanList);
                                    }
                                    dIA_Setting.DIAWindows.get(new XYData(precursormz - Loffset, precursormz + Roffset)).add(scanno);
                                    substr = substr.substring(substr.indexOf("</precursorMz>") + 14);
                                }
                            } else if (datatype == SpectralDataType.DataType.DIA_F_Window || datatype == SpectralDataType.DataType.pSMART || datatype == SpectralDataType.DataType.WiSIM) {
                                int stopidx = temp.indexOf("</precursorMz>");
                                if (stopidx == -1) {
                                    LogManager.getRootLogger().error("Parsing </precursorMz> failed. scan number :" + scanno);                                    
                                    System.exit(3);
                                }
                                int startidx = 0;
                                for (int i = stopidx; i > 0; i--) {
                                    if (temp.charAt(i) == '>') {
                                        startidx = i + 1;
                                        break;
                                    }
                                }
                                float precursormz = Float.parseFloat(temp.substring(startidx, stopidx));
                                //By default, assuming it's 5600 data, 
                                //and assume the precursor m/z is at 0.25 * window size Da to the lower bound of isolation window
                                float Loffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.2f;
                                float Roffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.8f;
                                
                                //If the scan contains "windowWideness", then it is a Thermo data, overwrite the isolation window ranges
                                if (temp.contains("windowWideness=\"")) {
                                    startidx = temp.indexOf("windowWideness=\"") + 16;
                                    stopidx = startidx + 1;
                                    for (int i = startidx + 1; i < temp.length(); i++) {
                                        if (temp.charAt(i) == '\"') {
                                            stopidx = i;
                                            break;
                                        }
                                    }
                                    float windowwideness = Float.parseFloat(temp.substring(startidx, stopidx));
                                     //Again assume the precursor m/z is at the center of isolation window, because it is a Thermo data
                                    Loffset = windowwideness / 2f;
                                    Roffset = windowwideness / 2f;
                                }

                                if (!dIA_Setting.DIAWindows.containsKey(new XYData(precursormz - Loffset, precursormz + Roffset))) {
                                    ArrayList<Integer> scanList = new ArrayList<>();
                                    dIA_Setting.DIAWindows.put(new XYData(precursormz - Loffset, precursormz + Roffset), scanList);
                                }
                                dIA_Setting.DIAWindows.get(new XYData(precursormz - Loffset, precursormz + Roffset)).add(scanno);
                            } else if (datatype == SpectralDataType.DataType.DIA_V_Window) {
                                //if the DIA data is variable window size setting, then use the pre-defined setting
                                int stopidx = temp.indexOf("</precursorMz>");
                                int startidx = 0;
                                for (int i = stopidx; i > 0; i--) {
                                    if (temp.charAt(i) == '>') {
                                        startidx = i + 1;
                                        break;
                                    }
                                }
                                float precursormz = Float.parseFloat(temp.substring(startidx, stopidx));
                                for (XYData window : dIA_Setting.DIAWindows.keySet()) {
                                    if (window.getX() <= precursormz && window.getY() >= precursormz) {
                                        dIA_Setting.DIAWindows.get(window).add(scanno);
                                        break;
                                    }
                                }
                            } else if (datatype == SpectralDataType.DataType.MSe) {
                                float mzlowF = 0f;
                                float mzhighF = 10000f;
                                if (!dIA_Setting.DIAWindows.containsKey(new XYData(mzlowF, mzhighF))) {
                                    ArrayList<Integer> scanList = new ArrayList<>();
                                    dIA_Setting.DIAWindows.put(new XYData(mzlowF, mzhighF), scanList);
                                }
                                dIA_Setting.DIAWindows.get(new XYData(mzlowF, mzhighF)).add(scanno);
                            }
                        }
                    }
                } else {
                    LogManager.getRootLogger().error("index of mzXML error");
                    System.exit(1);
                }
                ElutionTimeToScanNoMap.put(rt, scanno);
                ScanToElutionTime.put(scanno, rt);
                MsLevelList.put(scanno, mslevel);
            }
            ScanToElutionTime.compact();
            fileHandler.close();
        }
    }
     
     //Get all the DIA MS2 scans according to a isolation window range
    @Override
    public ScanCollection GetScanDIAMS2(XYData DIAWindow, boolean IncludePeak, float startTime, float endTime) {
        if (dIA_Setting == null) {
            LogManager.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
        ScanCollection swathScanCollection = new ScanCollection(parameter.Resolution);
        List<MzXMLthreadUnit> ScanList = new ArrayList<>();

        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
//        ArrayList<Integer> IncludedScans=new ArrayList<>();
        final BitSet IncludedScans=new BitSet();
        for(int scannum :dIA_Setting.DIAWindows.get(DIAWindow)){
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                IncludedScans.set(scannum, true);
            }
        }
        ScanList=ParseScans(IncludedScans);        
        for (MzXMLthreadUnit result : ScanList) {
            swathScanCollection.AddScan(result.scan);
            swathScanCollection.ElutionTimeToScanNoMap.put(result.scan.RetentionTime, result.scan.ScanNum);
        }        
        ScanList.clear();
        ScanList = null;
        return swathScanCollection;
    }
        
     //Get all the DIA MS1 scans according to MS1 m/z range, this was only for WiSIM data
    @Override
    public ScanCollection GetScanCollectionMS1Window(XYData MS1Window, boolean IncludePeak, float startTime, float endTime)  {
        if (dIA_Setting == null) {
            LogManager.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
        ScanCollection MS1WindowScanCollection = new ScanCollection(parameter.Resolution);
       
        List<MzXMLthreadUnit> ScanList = null;

        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
//        ArrayList<Integer> IncludedScans=new ArrayList<>();
        final BitSet IncludedScans=new BitSet();
        for(int scannum : dIA_Setting.MS1Windows.get(MS1Window)){
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                IncludedScans.set(scannum, true);
            }
        }
        
        ScanList=ParseScans(IncludedScans);
        
        for (MzXMLthreadUnit result : ScanList) {
            MS1WindowScanCollection.AddScan(result.scan);
            MS1WindowScanCollection.ElutionTimeToScanNoMap.put(result.scan.RetentionTime, result.scan.ScanNum);
        }
        ScanList.clear();
        ScanList = null;
        
        return MS1WindowScanCollection;
    }
    static private int step = -1;
    //Parse scans given a list of scan numbers
//    private List<MzXMLthreadUnit>  ParseScans(ArrayList<Integer> IncludedScans){
    private List<MzXMLthreadUnit> ParseScans(final BitSet IncludedScans){
         List<MzXMLthreadUnit> ScanList=new ArrayList<>();
         ArrayList<ForkJoinTask<?>> futures = new ArrayList<>();
        final ForkJoinPool fjp= new ForkJoinPool(NoCPUs);
        Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();        
        Entry<Integer, Long> ent = iter.next();
        long currentIdx = ent.getValue();
        int nextScanNo = ent.getKey();
        final RandomAccessFile fileHandler;
        try{fileHandler = new RandomAccessFile(filename, "r");}
        catch(FileNotFoundException e){throw new RuntimeException(e);}
        byte[] buffer = new byte[1<<10];
        if(step==-1)
            step=fjp.getParallelism()*32;
        while (iter.hasNext()) {
            ent = iter.next();
            long startposition = currentIdx;
            long nexposition = ent.getValue();
            int currentScanNo = nextScanNo;
            nextScanNo = ent.getKey();
            currentIdx = nexposition;

            if (IncludedScans.get(currentScanNo)) {
                try {
                    final int bufsize =  (int) (nexposition - startposition);
                    if(buffer.length<bufsize)
                        buffer = new byte[Math.max(bufsize,buffer.length<<1)];
//                    byte[] buffer = new byte[bufsize];
//                    RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
                    fileHandler.seek(startposition);
                    fileHandler.read(buffer, 0, bufsize);
//                    fileHandler.close();
//                    String xmltext = new String(buffer);
                    String xmltext = new String(buffer,0,bufsize,StandardCharsets.ISO_8859_1);
                    if (ent.getKey() == Integer.MAX_VALUE) {
                        xmltext = xmltext.replaceAll("</msRun>", "");
                    }
                    boolean ReadPeak = true;
                    final MzXMLthreadUnit unit = new MzXMLthreadUnit(xmltext, parameter, datatype, ReadPeak);
                    futures.add(fjp.submit(unit));
                    ScanList.add(unit);

                    if ((ScanList.size() % step) == 0) {
                        futures.get(futures.size()-step).get();
                        if (iter.hasNext() && fjp.getActiveThreadCount() < fjp.getParallelism()) {
                            step *= 2;
//                            System.out.println("MzXMLthreadUnit: fjp.getActiveThreadCount()\t" + fjp.getActiveThreadCount()+"\t"+step);
                        }
                    }
                } catch (Exception ex) {
                    LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                }
            }
        }
        try {fileHandler.close();}
        catch (IOException ex) {throw new RuntimeException(ex);}
        fjp.shutdown();
        try {fjp.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);}
        catch (InterruptedException ex) {throw new RuntimeException(ex);}
//        for (MzXMLthreadUnit unit : ScanList) {
//            executorPool.execute(unit);
//        }
//        executorPool.shutdown();
//
//        try {
//            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//        } catch (InterruptedException e) {
//            LogManager.getRootLogger().info("interrupted..");
//        }
        return ScanList;
    }
    
    @Override
    public ScanCollection GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak, float startTime, float endTime) {
        ScanCollection scanCollection = InitializeScanCollection();
        LogManager.getRootLogger().debug("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB (" + NoCPUs + " threads)");

//        ArrayList<Integer> IncludedScans = new ArrayList<>();
        final IntArrayList IncludedScans = new IntArrayList();
        
        for(int ScanNum : MsLevelList.keySet()){
            if(MsLevelList.get(ScanNum)==1 && MS1Included){
                IncludedScans.add(ScanNum);
            }
            if(MsLevelList.get(ScanNum)==2 && MS2Included){
                IncludedScans.add(ScanNum);
            }
        }
         
        List<MzXMLthreadUnit> ScanList = null;

        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
        
//        ArrayList<Integer> temp=new ArrayList<>();
        final BitSet temp = new BitSet();
//        for(int scannum : IncludedScans){
        for(int i=0; i<IncludedScans.size(); ++i){
            final int scannum = IncludedScans.get(i);
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                temp.set(scannum, true);
            }
        }
        
        ScanList=ParseScans(temp);
        
        for (MzXMLthreadUnit result : ScanList) {
            scanCollection.AddScan(result.scan);
        }
        ScanList.clear();
        ScanList = null;
        
        System.gc();
        LogManager.getRootLogger().debug("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");
        return scanCollection;
    }
 
}
