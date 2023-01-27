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
package MSUmpire.PSMDataStructure;

import com.compomics.util.experiment.biology.Enzyme;
import com.compomics.util.experiment.biology.EnzymeFactory;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

/**
 * Enzyme manager from compomics library
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class EnzymeManager {

    private static EnzymeManager enzymeManager;
    EnzymeFactory enzymeFactory = null;

    public EnzymeManager() throws IOException {
        if (enzymeFactory == null) {
            enzymeFactory = EnzymeFactory.getInstance();
            InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/enzymes.xml");
            String tmpfile = "enzymes.xml";
            InputStreamToFile convert = new InputStreamToFile();

            File enzymeFile = convert.GetFile(is, tmpfile);
            try {
                enzymeFactory.importEnzymes(enzymeFile);
            } catch (Exception e) {
                throw new RuntimeException("Could not import enzymes");
            }
        }
    }
    
    public static EnzymeManager GetInstance() throws IOException {
        if (enzymeManager == null) {
            enzymeManager = new EnzymeManager();
        }
        return enzymeManager;
    }

    public Enzyme GetTrypsin() {
        return enzymeFactory.getEnzyme("Trypsin");
    }
        
    public Enzyme GetTrypsinNoP() {
        return enzymeFactory.getEnzyme("Trypsin, no P rule");
    }
    
    public Enzyme GetSemiTryptic(){
        return enzymeFactory.getEnzyme("Semi-Tryptic");        
    }
}
