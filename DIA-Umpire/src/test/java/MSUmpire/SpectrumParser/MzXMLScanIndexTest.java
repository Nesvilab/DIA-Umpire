package MSUmpire.SpectrumParser;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Base64;

import static org.junit.Assert.assertEquals;

public class MzXMLScanIndexTest {
    @Rule public TemporaryFolder temporaryFolder = new TemporaryFolder();

    @Test
    public void readsContiguousScanNumbers() throws Exception {
        assertWindowRouting(new int[]{1, 2, 3});
    }

    @Test
    public void readsScanNumbersWithGap() throws Exception {
        assertWindowRouting(new int[]{1, 3, 4});
    }

    private void assertWindowRouting(int[] scanNumbers) throws Exception {
        Path file = temporaryFolder.newFile("scans-" + scanNumbers[1] + ".mzXML").toPath();
        StringBuilder xml = new StringBuilder("<?xml version=\"1.0\"?>\n<mzXML>\n<msRun>\n");
        for (int i = 0; i < scanNumbers.length; ++i) {
            String peaks = Base64.getEncoder().encodeToString(ByteBuffer.allocate(16)
                    .putFloat(200 + i).putFloat(1000).putFloat(300 + i).putFloat(500).array());
            xml.append("<scan num=\"").append(scanNumbers[i]).append("\" msLevel=\"")
                    .append(i == 0 ? 1 : 2)
                    .append("\" peaksCount=\"2\" centroided=\"1\" retentionTime=\"PT")
                    .append(i + 1).append("S\">\n");
            if (i > 0) {
                xml.append("<precursorMz precursorIntensity=\"1\" windowWideness=\"2\">")
                        .append(i == 1 ? 500 : 700).append("</precursorMz>\n");
            }
            xml.append("<peaks precision=\"32\" byteOrder=\"network\" pairOrder=\"m/z-int\">")
                    .append(peaks).append("</peaks>\n</scan>\n");
        }
        Files.write(file, xml.append("</msRun>\n</mzXML>\n").toString()
                .getBytes(StandardCharsets.US_ASCII));

        InstrumentParameter parameters = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600);
        parameters.EstimateBG = false;
        parameters.MinMSMSIntensity = 0;
        DIA_Setting setting = new DIA_Setting();
        mzXMLParser parser = new mzXMLParser(file.toString(), parameters,
                SpectralDataType.DataType.DIA_F_Window, setting, 1);

        assertEquals(2, setting.DIAWindows.size());
        for (XYData window : setting.DIAWindows.keySet()) {
            int expectedScan = window.getX() < 600 ? scanNumbers[1] : scanNumbers[2];
            ScanCollection scans = parser.GetScanDIAMS2(window, true, 0, 999);
            assertEquals(1, scans.ScanHashMap.size());
            assertEquals(expectedScan, scans.ScanHashMap.values().iterator().next().ScanNum);
        }
    }
}
