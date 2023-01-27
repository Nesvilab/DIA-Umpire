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

import MSUmpire.BaseDataStructure.ScanData;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.StringReader;
import java.nio.ByteBuffer;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.codec.binary.Base64;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

//import uk.ac.ebi.pride.jaxb.utils.BinaryDataUtils;
//import uk.ac.ebi.pride.jaxb.utils.CvTermReference;
/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class mzXMLReadUnit {

    String XMLtext;

    public mzXMLReadUnit(String XMLtext) {
        this.XMLtext = XMLtext;
    }
    private static synchronized DocumentBuilder get_docBuilder() {
        try{return DocumentBuilderFactory.newInstance().newDocumentBuilder();}
        catch(ParserConfigurationException ex){throw new RuntimeException(ex);}
    }

    private static final ThreadLocal<DocumentBuilder> tls= new ThreadLocal<DocumentBuilder>(){
        @Override
        protected DocumentBuilder initialValue(){
            return get_docBuilder();
        }
    };

    public ScanData Parse() throws ParserConfigurationException, SAXException, IOException, DataFormatException {
        if (XMLtext.replaceFirst("</scan>", "").contains("</scan>")) {
            XMLtext = XMLtext.replaceFirst("</scan>", "");
        }
        if (!XMLtext.contains("</scan>")) {
            XMLtext += "</scan>";
        }
        ScanData scan = new ScanData();
        final DocumentBuilder docBuilder = tls.get();
        docBuilder.reset();
        InputSource input = new InputSource(new StringReader(XMLtext));
        Document doc = null;
        try {
            doc = docBuilder.parse(input);
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            LogManager.getRootLogger().error(XMLtext);
        }
        Node root = doc.getFirstChild();
        for (int i = 0; i < root.getAttributes().getLength(); i++) {
            switch (root.getAttributes().item(i).getNodeName()) {
                case ("num"):
                    scan.ScanNum = Integer.parseInt(root.getAttributes().item(i).getNodeValue());
                    break;
                case ("centroided"): {
                    if ("1".equals(root.getAttributes().item(i).getNodeValue())) {
                        scan.centroided = true;
                    } else {
                        scan.centroided = false;
                    }
                    break;
                }
                case ("msLevel"):
                    scan.MsLevel = Integer.parseInt(root.getAttributes().item(i).getNodeValue());
                    break;
                case ("scanType"):
                    scan.scanType = root.getAttributes().item(i).getNodeValue();
                    break;
                case ("peaksCount"):
                    scan.PeaksCountString = Integer.parseInt(root.getAttributes().item(i).getNodeValue());
                    break;
                case ("retentionTime"):
                    scan.RetentionTime = Float.parseFloat(root.getAttributes().item(i).getNodeValue().substring(2, root.getAttributes().item(i).getNodeValue().indexOf("S"))) / 60f;
                    break;
                case ("lowMz"): {
                    String value = root.getAttributes().item(i).getNodeValue();
                    if ("inf".equals(value)) {
                        value = String.valueOf(Float.MIN_VALUE);
                    }
                    scan.StartMz = Float.parseFloat(value);
                    break;
                }
                case ("highMz"): {
                    String value = root.getAttributes().item(i).getNodeValue();
                    if ("inf".equals(value)) {
                        value = String.valueOf(Float.MAX_VALUE);
                    }
                    scan.EndMz = Float.parseFloat(value);
                    break;
                }
                case ("startMz"):
                    scan.StartMz = Float.parseFloat(root.getAttributes().item(i).getNodeValue());
                    break;
                case ("endMz"): {
                    String value = root.getAttributes().item(i).getNodeValue();
                    if ("inf".equals(value)) {
                        value = String.valueOf(Float.MAX_VALUE);
                    }
                    scan.EndMz = Float.parseFloat(value);
                    break;
                }
                case ("basePeakMz"): {
                    String value = root.getAttributes().item(i).getNodeValue();
                    if ("inf".equals(value)) {
                        value = String.valueOf(Float.MAX_VALUE);
                    }
                    if (!"null".contentEquals(value))
                        scan.BasePeakMz = Float.parseFloat(value);
                    break;
                }
                case ("basePeakIntensity"):
                    final String value = root.getAttributes().item(i).getNodeValue();
                    if (!"null".contentEquals(value))
                        scan.BasePeakIntensity = Float.parseFloat(value);
                    break;
                case ("totIonCurrent"):
                    scan.SetTotIonCurrent(Float.parseFloat(root.getAttributes().item(i).getNodeValue()));
                    break;
            }
        }
        for (int i = 0; i < root.getChildNodes().getLength(); i++) {
            Node childNode = root.getChildNodes().item(i);
            switch (childNode.getNodeName()) {
                case ("precursorMz"): {
                    scan.PrecursorMz = Float.parseFloat(childNode.getTextContent());
                    for (int j = 0; j < childNode.getAttributes().getLength(); j++) {
                        switch (childNode.getAttributes().item(j).getNodeName()) {
                            case ("precursorScanNum"):
                                scan.precursorScanNum = Integer.parseInt(childNode.getAttributes().item(j).getNodeValue());
                                break;
                            case ("precursorIntensity"):
                                scan.PrecursorIntensity = Float.parseFloat(childNode.getAttributes().item(j).getNodeValue());
                                break;
                            case ("precursorCharge"):
                                scan.PrecursorCharge = Integer.parseInt(childNode.getAttributes().item(j).getNodeValue());
                                break;
                            case ("activationMethod"):
                                scan.ActivationMethod = childNode.getAttributes().item(j).getNodeValue();
                                break;
                            case ("windowWideness"):
                                scan.windowWideness = Float.parseFloat(childNode.getAttributes().item(j).getNodeValue());
                                break;
                        }
                    }
                    break;
                }
                case ("peaks"): {
                    for (int j = 0; j < childNode.getAttributes().getLength(); j++) {
                        switch (childNode.getAttributes().item(j).getNodeName()) {
                            case ("compressionType"):
                                scan.compressionType = childNode.getAttributes().item(j).getNodeValue();
                                break;
                            case ("precision"):
                                scan.precision = Integer.parseInt(childNode.getAttributes().item(j).getNodeValue());
                                break;
                        }
                    }
                    ParsePeakString(scan, childNode.getTextContent());
                    break;
                }
            }
            childNode = null;
        }
        if ("calibration".equals(scan.scanType)) {
            scan.MsLevel = -1;
        }
        XMLtext = null;
        scan.Data.Finalize();
        return scan;
    }

    public byte[] ZlibUncompressBuffer(byte[] compressed) throws IOException, DataFormatException {

        Inflater decompressor = new Inflater();
        decompressor.setInput(compressed);

        ByteArrayOutputStream bos = null;
        try {

            bos = new ByteArrayOutputStream(compressed.length);

            // Decompress the data
            byte[] buf = new byte[decompressor.getRemaining() * 2];
            while (decompressor.getRemaining() > 0) {
                int count = decompressor.inflate(buf);
                bos.write(buf, 0, count);
            }

        } finally {
            try {
                bos.close();
            } catch (Exception nope) { /* This exception doesn't matter */ }
        }
        decompressor.end();
        compressed = null;
        decompressor = null;
        byte[] result = bos.toByteArray();
        bos = null;
        return result;
    }

    private void ParsePeakString(ScanData scan, String peakString) throws IOException, DataFormatException {
        int offset;

        peakString = peakString.replaceAll("\n", "");
        byte[] decoded = Base64.decodeBase64(peakString.getBytes());

        if ("zlib".equals(scan.compressionType)) {
            decoded = ZlibUncompressBuffer(decoded);
        }
        switch (scan.precision) {
            case (32): {
                offset = 0;
                for (int i = 0; i < scan.PeaksCountString; i++) {
                    byte[] mz = new byte[]{decoded[offset], decoded[offset + 1], decoded[offset + 2], decoded[offset + 3]};
                    byte[] intensity = new byte[]{decoded[offset + 4], decoded[offset + 5], decoded[offset + 6], decoded[offset + 7]};
                    ByteBuffer mzBuffer = ByteBuffer.wrap(mz);
                    ByteBuffer intBuffer = ByteBuffer.wrap(intensity);
                    float intensityfloat = intBuffer.getFloat();
                    float mzfloat = mzBuffer.getFloat();
                    if (intensityfloat > 0f) {
                        scan.AddPoint(mzfloat, intensityfloat);
                    }
                    mz = null;
                    intensity = null;
                    mzBuffer.clear();
                    intBuffer.clear();
                    mzBuffer = null;
                    intBuffer = null;
                    offset += 8;
                }
                break;
            }
            case (64): {
                offset = 0;
                for (int i = 0; i < scan.PeaksCountString; i++) {
                    byte[] mz = new byte[]{decoded[offset], decoded[offset + 1], decoded[offset + 2], decoded[offset + 3], decoded[offset + 4], decoded[offset + 5], decoded[offset + 6], decoded[offset + 7]};
                    byte[] intensity = new byte[]{decoded[offset + 8], decoded[offset + 9], decoded[offset + 10], decoded[offset + 11], decoded[offset + 12], decoded[offset + 13], decoded[offset + 14], decoded[offset + 15]};
                    ByteBuffer mzBuffer = ByteBuffer.wrap(mz);
                    ByteBuffer intBuffer = ByteBuffer.wrap(intensity);
                    float intensityfloat = (float) intBuffer.getDouble();
                    float mzfloat = (float) mzBuffer.getDouble();
                    if (intensityfloat > 0f) {
                        scan.AddPoint(mzfloat, intensityfloat);
                    }
                    mz = null;
                    intensity = null;
                    mzBuffer.clear();
                    intBuffer.clear();
                    mzBuffer = null;
                    intBuffer = null;
                    offset += 16;
                }
                break;
            }
        }
        peakString = null;
        decoded = null;
    }
}
