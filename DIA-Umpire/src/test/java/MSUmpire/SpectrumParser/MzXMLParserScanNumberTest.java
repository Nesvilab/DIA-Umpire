package MSUmpire.SpectrumParser;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.Deflater;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

/**
 * mzXMLParser on DIA files written in the to_mzXML layout. The run has two isolation windows, and one MS2 spectrum
 * (num 5) is absent, as when to_mzXML drops a spectrum without peaks.
 */
public class MzXMLParserScanNumberTest {

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    private static final float ISOLATION_WIDTH = 150f;
    private static final float LOW_TARGET = 425f;   // isolates 350-500
    private static final float HIGH_TARGET = 575f;  // isolates 500-650
    private static final float MS1 = 0f;

    private static final int[] NUMS = {1, 2, 3, 4, 6, 7, 8, 9};
    private static final float[] TARGETS = {MS1, LOW_TARGET, HIGH_TARGET, MS1, HIGH_TARGET, MS1, LOW_TARGET, HIGH_TARGET};

    @Test
    public void scanIndexKeysAreTheNumAttributes() throws Exception {
        final mzXMLParser parser = Open(WriteRun("keys", NUMS, TARGETS), new DIA_Setting());

        final List<Integer> keys = new ArrayList<>(parser.ScanIndex.keySet());
        assertEquals(Integer.MAX_VALUE, (int) keys.remove(keys.size() - 1));
        assertEquals(Arrays.asList(1, 2, 3, 4, 6, 7, 8, 9), keys);
    }

    @Test
    public void eachIsolationWindowGetsItsOwnScans() throws Exception {
        final DIA_Setting dia = new DIA_Setting();
        final mzXMLParser parser = Open(WriteRun("windows", NUMS, TARGETS), dia);

        final TreeMap<Float, TreeSet<Integer>> numsByTarget = new TreeMap<>();
        for (int i = 0; i < NUMS.length; ++i) {
            if (TARGETS[i] != MS1) {
                numsByTarget.computeIfAbsent(TARGETS[i], t -> new TreeSet<>()).add(NUMS[i]);
            }
        }
        final TreeSet<Float> targetsSeen = new TreeSet<>();
        for (Map.Entry<XYData, ArrayList<Integer>> isolation : dia.DIAWindows.entrySet()) {
            final float target = (isolation.getKey().getX() + isolation.getKey().getY()) / 2f;
            targetsSeen.add(target);
            final TreeSet<Integer> listed = new TreeSet<>(isolation.getValue());
            assertEquals(numsByTarget.get(target), listed);

            final ScanCollection read = parser.GetScanDIAMS2(isolation.getKey(), true, 0f, 1000f);
            assertEquals("spectra read for target " + target, listed, new TreeSet<>(read.ScanHashMap.keySet()));
            for (ScanData spectrum : read.ScanHashMap.values()) {
                assertEquals(2, spectrum.MsLevel);
                assertEquals(target, spectrum.PrecursorMz, 1e-3f);
                // The only peak encodes the num, so this checks the bytes read, not only the scan tag.
                assertEquals(1, spectrum.PointCount());
                assertEquals(PeakMz(spectrum.ScanNum), spectrum.Data.get(0).getX(), 1e-3f);
                assertEquals(PeakIntensity(spectrum.ScanNum), spectrum.Data.get(0).getY(), 1e-1f);
            }
        }
        assertEquals(numsByTarget.keySet(), targetsSeen);
    }

    @Test
    public void stopsWhenScanNumbersDoNotIncrease() throws Exception {
        final File run = WriteRun("decreasing", new int[]{1, 3, 2}, new float[]{MS1, LOW_TARGET, HIGH_TARGET});

        final IOException error = assertThrows(IOException.class, () -> Open(run, new DIA_Setting()));
        assertTrue(error.getMessage(), error.getMessage().contains("scan number 2"));
    }

    private static float PeakMz(final int num) {
        return 150.5f + num;
    }

    private static float PeakIntensity(final int num) {
        return 10f + 50f * num;
    }

    private static mzXMLParser Open(final File mzXML, final DIA_Setting dia) throws Exception {
        final InstrumentParameter parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.Orbitrap);
        parameter.EstimateBG = false;
        parameter.Denoise = false;
        return new mzXMLParser(mzXML.getPath(), parameter, SpectralDataType.DataType.DIA_F_Window, dia, 2);
    }

    // One spectrum per num, each with a single peak, in the text layout of mzXMLParser.to_mzXML.
    private File WriteRun(final String name, final int[] nums, final float[] targets) throws IOException {
        final StringBuilder out = new StringBuilder("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<mzXML>\n  <msRun>\n");
        for (int i = 0; i < nums.length; ++i) {
            final int num = nums[i];
            out.append(String.format(Locale.ROOT, "    <scan num=\"%d\"\n"
                    + "          scanType=\"Full\"\n"
                    + "          centroided=\"1\"\n"
                    + "          msLevel=\"%d\"\n"
                    + "          peaksCount=\"1\"\n"
                    + "          polarity=\"+\"\n"
                    + "          retentionTime=\"PT%sS\"\n"
                    + "          lowMz=\"%f\"\n"
                    + "          highMz=\"%f\"\n"
                    + "          basePeakMz=\"%f\"\n"
                    + "          basePeakIntensity=\"%f\"\n"
                    + "          totIonCurrent=\"%f\">\n",
                    num, targets[i] == MS1 ? 1 : 2, Float.toString(3f * i + 1f),
                    PeakMz(num), PeakMz(num), PeakMz(num), PeakIntensity(num), PeakIntensity(num)));
            if (targets[i] != MS1) {
                out.append(String.format(Locale.ROOT, "      <precursorMz precursorIntensity=\"0.000000\" activationMethod=\"HCD\" windowWideness=\"%f\">%f</precursorMz>\n",
                        ISOLATION_WIDTH, targets[i]));
            }
            final byte[] raw = ByteBuffer.allocate(2 * Float.BYTES).putFloat(PeakMz(num)).putFloat(PeakIntensity(num)).array();
            final Deflater zlib = new Deflater();
            zlib.setInput(raw);
            zlib.finish();
            final byte[] packed = new byte[raw.length + 64];
            final int packedLength = zlib.deflate(packed);
            zlib.end();
            out.append(String.format(Locale.ROOT, "      <peaks compressionType=\"zlib\"\n"
                    + "             compressedLen=\"%d\"\n"
                    + "             precision=\"32\"\n"
                    + "             byteOrder=\"network\"\n"
                    + "             contentType=\"m/z-int\">%s</peaks>\n"
                    + "    </scan>\n",
                    packedLength, Base64.getEncoder().encodeToString(Arrays.copyOf(packed, packedLength))));
        }
        out.append("  </msRun>\n</mzXML>");
        final File file = new File(tmp.newFolder(name), name + ".mzXML");
        Files.write(file.toPath(), out.toString().getBytes(StandardCharsets.ISO_8859_1));
        return file;
    }
}
