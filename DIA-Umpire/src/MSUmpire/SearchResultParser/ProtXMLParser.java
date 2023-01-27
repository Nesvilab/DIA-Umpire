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
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import org.apache.logging.log4j.LogManager;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ProtXMLParser {

    public LCMSID SingleLCMSID;
    public String FileName;
    public float threshold = 0f;

    public ProtXMLParser(LCMSID singleLCMSID, String FileName, float threshold) throws IOException, ClassNotFoundException, InterruptedException {
        this.SingleLCMSID = singleLCMSID;
        this.FileName = FileName;
        this.threshold = threshold;
        LogManager.getRootLogger().info("Parsing protXML: " + FileName + ".....");
        ParseSAX();
        SingleLCMSID.DetermineAssignIonListByProtPepSeq();
        SingleLCMSID.FixProteinWithDecoyHead();
        singleLCMSID.UpdateDecoyMaxIniProb();
        SingleLCMSID.SetGroupProbForNonDecoyGroupHead();
        //System.out.print("done\n");
    }


   private void ParseSAX() throws IOException {
        File fXmlFile = new File(FileName);
        if (!fXmlFile.exists()) {
            LogManager.getRootLogger().info("File :" + FileName + " cannot be found\n");
            return;
        }
        FileInputStream inputStream = new FileInputStream(FileName);
        SolnaParser parser = new SolnaParser();
        ProtXMLParseHandler handler = new ProtXMLParseHandler(SingleLCMSID,threshold);
        parser.addHandler("/protein_summary/protein_group", handler);
        parser.parse(inputStream);        
    }    
}
