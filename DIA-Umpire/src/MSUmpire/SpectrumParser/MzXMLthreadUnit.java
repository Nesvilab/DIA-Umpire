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
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.xml.sax.SAXException;

/**
 * Thread unit for parsing one scan in mzXML file
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MzXMLthreadUnit implements Runnable {

    public ScanData scan;
    private String XMLtext;
    private InstrumentParameter parameter;
    boolean ReadPeak = true;
    SpectralDataType.DataType dataType = SpectralDataType.DataType.DDA;

    public MzXMLthreadUnit(String XMLtext, InstrumentParameter parameter, SpectralDataType.DataType dataType,boolean ReadPeak) {
        this.XMLtext = XMLtext;
        this.parameter = parameter;
        this.ReadPeak = ReadPeak;
        this.dataType = dataType;
    }

    public MzXMLthreadUnit(String XMLtext, InstrumentParameter parameter, SpectralDataType.DataType dataType) {
        this.XMLtext = XMLtext;
        this.parameter = parameter;
        this.dataType = dataType;
    }

    private void Read() throws FileNotFoundException, IOException, ParserConfigurationException, SAXException, DataFormatException {
        mzXMLReadUnit read = new mzXMLReadUnit(this.XMLtext);
        this.scan = read.Parse();
        this.XMLtext = null;
        read = null;
    }

    @Override
    public void run() {
        try {
            Read();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
        
       scan.Preprocessing(parameter);
       XMLtext=null;
       parameter=null;
    }
}
