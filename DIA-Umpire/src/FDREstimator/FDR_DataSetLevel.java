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
package FDREstimator;

import MSUmpire.BaseDataStructure.DBSearchParam;
import MSUmpire.DIA.DIAPack;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import java.io.IOException;
import java.util.ArrayList;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FDR_DataSetLevel {

    public LCMSID combineID = null;

    public void GeneratePepIonList(ArrayList<DIAPack> DIAFileList, DBSearchParam param, String combineIDPath) throws IOException, ParserConfigurationException, SAXException, ClassNotFoundException, InterruptedException {

        for (DIAPack diafile : DIAFileList) {
            diafile.ParsePepXML(param,null);
        }

        //Estimate peptide level PepFDR in whole dataset
        combineID = new LCMSID(combineIDPath, param.DecoyPrefix, param.FastaPath);
        for (DIAPack Diafile : DIAFileList) {
            LCMSID lcms = Diafile.IDsummary;
            for (PepIonID pepIonID : lcms.GetPepIonList().values()) {
                if (!combineID.GetPepIonList().containsKey(pepIonID.GetKey())) {
                    PepIonID newpep = pepIonID.ClonePepIonID();
                    if (pepIonID.IsDecoy(param.DecoyPrefix)) {
                        newpep.IsDecoy = 1;
                    } else {
                        newpep.IsDecoy = 0;
                    }
                    combineID.AddPeptideID(newpep);
                }
                if (combineID.GetPepIonList().get(pepIonID.GetKey()).MaxProbability < pepIonID.MaxProbability) {
                    combineID.GetPepIonList().get(pepIonID.GetKey()).MaxProbability = pepIonID.MaxProbability;
                }
            }
        }
        combineID.DecoyTag = param.DecoyPrefix;
        combineID.FDR = param.PepFDR;
        combineID.FindPepProbThresholdByFDR();
        combineID.RemoveDecoyPep();
        combineID.RemoveLowProbPep();
    }
}
