package MSUmpire.SpectrumParser;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import org.junit.Test;

public class InferMzXMLOffsetsTest {

    private static Infer_mzXML_offsets.ScanOffsets Infer(final String xml) throws IOException {
        return Infer_mzXML_offsets.infer_mzXML_scan_offsets(new ByteArrayInputStream(xml.getBytes(StandardCharsets.ISO_8859_1)));
    }

    @Test
    public void returnsTheNumAndOffsetOfEachScanTag() throws IOException {
        final String xml = "<mzXML>\n  <msRun>\n"
                + "    <scan num=\"10\"\n          msLevel=\"1\">\n      <peaks>AAAA</peaks>\n    </scan>\n"
                + "    <scan num=\"12\"\n          msLevel=\"2\">\n      <peaks>AAAA</peaks>\n    </scan>\n"
                + "    <scan num=\"17\" msLevel=\"2\"><peaks>AAAA</peaks></scan>\n"
                + "  </msRun>\n</mzXML>";
        final Infer_mzXML_offsets.ScanOffsets scans = Infer(xml);

        assertArrayEquals(new int[]{10, 12, 17}, scans.scanNums.toArray());
        assertEquals(4, scans.offsets.size());
        for (int i = 0; i < scans.scanNums.size(); ++i) {
            final int offset = (int) scans.offsets.get(i);
            assertEquals("<scan num=\"" + scans.scanNums.get(i) + "\"", xml.substring(offset, xml.indexOf('"', xml.indexOf('"', offset) + 1) + 1));
        }
        final int end = (int) scans.offsets.getLast();
        assertEquals("</scan>", xml.substring(end - "</scan>".length(), end));
    }

    @Test
    public void findsSeveralScanTagsOnOneLine() throws IOException {
        final String xml = "<msRun><scan num=\"7\" msLevel=\"1\"></scan><scan num=\"8\" msLevel=\"2\"></scan>\n<scan num=\"11\" msLevel=\"2\">\n</scan>\n</msRun>";
        final Infer_mzXML_offsets.ScanOffsets scans = Infer(xml);

        assertArrayEquals(new int[]{7, 8, 11}, scans.scanNums.toArray());
        assertArrayEquals(new long[]{xml.indexOf("<scan num=\"7\""), xml.indexOf("<scan num=\"8\""), xml.indexOf("<scan num=\"11\""),
            xml.lastIndexOf("</scan>") + "</scan>".length()}, scans.offsets.toArray());
    }

    @Test
    public void skipsElementsWhoseNamesStartWithScan() throws IOException {
        final String xml = "<msRun>\n"
                + "<scan num=\"3\" msLevel=\"1\">\n<scanOrigin parentFileID=\"x\" num=\"99\"/>\n<peaks>AAAA</peaks>\n</scan>\n"
                + "<scan num=\"4\" msLevel=\"1\"><peaks>AAAA</peaks>\n</scan>\n"
                + "</msRun>";
        final Infer_mzXML_offsets.ScanOffsets scans = Infer(xml);

        assertArrayEquals(new int[]{3, 4}, scans.scanNums.toArray());
        assertEquals(3, scans.offsets.size());
    }

    @Test
    public void ignoresAttributesThatEndInNum() throws IOException {
        final Infer_mzXML_offsets.ScanOffsets scans = Infer("<scan precursorScanNum=\"1\" num=\"2\" msLevel=\"2\">\n</scan>\n");

        assertArrayEquals(new int[]{2}, scans.scanNums.toArray());
    }

    @Test
    public void rejectsScanTagWithoutNum() {
        assertThrows(IOException.class, () -> Infer("<scan msLevel=\"1\">\n</scan>\n"));
    }

    @Test
    public void rejectsScanNumberReservedForEndOfIndex() {
        assertThrows(IOException.class, () -> Infer("<scan num=\"" + Integer.MAX_VALUE + "\">\n</scan>\n"));
    }
}
