/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
