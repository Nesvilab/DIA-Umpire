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
package MSUmpire.SearchResultParser;

import MSUmpire.PSMDataStructure.LCMSID;
import com.vseravno.solna.SolnaParser;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepXMLParser {

    private LCMSID singleLCMSID;
    public String FileName;
    public float threshold = 0f;
    public float StartRT = 0f;
    public float EndRT = 9999f;
    public boolean FilteredID = false;
    public boolean CorrectMassDiff=true;
    
    public PepXMLParser(LCMSID singleLCMSID, String FileName, float threshold, float StartRT, float EndRT) throws ParserConfigurationException, SAXException, IOException {
        this.singleLCMSID = singleLCMSID;
        this.FileName = FileName;
        this.threshold = threshold;
        this.StartRT = StartRT;
        this.EndRT = EndRT;
        LogManager.getRootLogger().info("Parsing pepXML: " + FileName + "....");
        ParseSAX();
    }

    
    public PepXMLParser(LCMSID singleLCMSID, String FileName, float threshold, boolean CorrectMassDiff) throws ParserConfigurationException, SAXException, IOException {
        this.singleLCMSID = singleLCMSID;
        this.CorrectMassDiff=CorrectMassDiff;
        this.FileName = FileName;
        this.threshold = threshold;
        LogManager.getRootLogger().info("Parsing pepXML: " + FileName + "....");
        try {
            ParseSAX();
        } catch (Exception e) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(e));
            LogManager.getRootLogger().info("Parsing pepXML: " + FileName + " failed. Trying to fix the file...");
            insert_msms_run_summary(new File(FileName));
            ParseSAX();
        }
        //System.out.print("done\n");
    }
    
    public PepXMLParser(LCMSID singleLCMSID, String FileName, float threshold) throws ParserConfigurationException, SAXException, IOException {
        this.singleLCMSID = singleLCMSID;
        this.FileName = FileName;
        this.threshold = threshold;
        LogManager.getRootLogger().info("Parsing pepXML: " + FileName + "....");
        try {
            ParseSAX();
        } catch (Exception e) {
            LogManager.getRootLogger().info("Parsing pepXML: " + FileName + " failed. Trying to fix the file...");
            insert_msms_run_summary(new File(FileName));
            ParseSAX();
        }
        //System.out.print("done\n");
    }

    private void ParseSAX() throws IOException {
        File fXmlFile = new File(FileName);
        if (!fXmlFile.exists()) {
            LogManager.getRootLogger().error("File :" + FileName + " cannot be found");
            return;
        }
        FileInputStream inputStream = new FileInputStream(FileName);
        SolnaParser parser = new SolnaParser();
        PepXMLParseHandler handler = new PepXMLParseHandler(singleLCMSID, StartRT, EndRT, threshold,CorrectMassDiff);
        //handler.FileBaseNameFilter=FilenameUtils.getBaseName(singleLCMSID.mzXMLFileName);
        parser.addHandler("/msms_pipeline_analysis/msms_run_summary/spectrum_query", handler);
        parser.addHandler("/msms_pipeline_analysis/msms_run_summary/search_summary", handler);
        parser.parse(inputStream);
    }

    public void insert_msms_run_summary(File inFile) throws IOException {
        // temp file
        File outFile = new File(FileName.replace("pep.xml", "_fixed.pep.xml"));
        FileName = FileName.replace("pep.xml", "_fixed.pep.xml");

        // input
        FileInputStream fis = new FileInputStream(inFile);
        BufferedReader in = new BufferedReader(new InputStreamReader(fis));

        // output         
        FileOutputStream fos = new FileOutputStream(outFile);
        PrintWriter out = new PrintWriter(fos);

        String thisLine = "";
        String lastLine = "";

        while ((thisLine = in.readLine()) != null) {
            if (!"".equals(lastLine)) {
                out.println(lastLine);
            }
            lastLine = thisLine;
        }
        if (!"".equals(lastLine)) {
            out.println("</msms_run_summary>");
            out.println(lastLine);
        }
        out.flush();
        out.close();
        in.close();
    }
}
