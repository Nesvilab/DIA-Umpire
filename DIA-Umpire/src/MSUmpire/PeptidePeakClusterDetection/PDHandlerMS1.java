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
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.LCMSPeakStructure.LCMSPeakBase;
import java.io.*;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;

/**
 * Peak detection processing class for MS1 peak
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PDHandlerMS1 extends PDHandlerBase {

    public PDHandlerMS1(LCMSPeakBase lcmspeak, int NoCPUs, float PPM) {
        this.NoCPUs = NoCPUs;
        this.PPM = PPM;
        this.LCMSPeakBase = lcmspeak;
        this.parameter = lcmspeak.parameter;
    }

    //Detect peak curve and isotope peak cluster given a list of ScanCollection
    public void DetectPeakClusters(ArrayList<ScanCollection> scanCollections) throws InterruptedException, ExecutionException, IOException {        
        FindAllMzTracePeakCurvesForScanCollections(scanCollections);
        scanCollections.clear();//not needed after this line, so clear it
        PeakCurveCorrClustering(new XYData(Float.NEGATIVE_INFINITY, Float.POSITIVE_INFINITY));
    }
}
