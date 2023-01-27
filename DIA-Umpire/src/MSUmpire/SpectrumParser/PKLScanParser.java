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
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PKLScanParser {

    public String filename;
    public ScanData scan;

    public PKLScanParser(String filename) throws IOException {
        this.filename = filename;
        Parse();
    }

    private void Parse() throws FileNotFoundException, IOException {
        //806.080993652344,8429.974609375,1
        //832.287536621094,7226.927734375,1
        //854.039978027344,6682.37646484375,1
        //861.061340332031,8370.4716796875,1
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String[] Values = null;
        scan = new ScanData();
        scan.MGFTitle = FilenameUtils.getBaseName(filename);
        while ((line = reader.readLine()) != null) {
            if ((Values = line.split(",")).length == 3) {
                scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
            }
        }
        reader.close();
    }
}
