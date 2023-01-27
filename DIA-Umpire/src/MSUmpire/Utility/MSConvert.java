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
package MSUmpire.Utility;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class MSConvert {

    public String msconvertpath="C:/inetpub/tpp-bin/msconvert";
    public String SpectrumPath;
    public MSConvert(String SpectrumPath){
        this.SpectrumPath=SpectrumPath;
    }
    public void Convert(){
        try {
            String[] msconvertcmd = {msconvertpath, "--mzXML", "--32", "-z", SpectrumPath, "-o", FilenameUtils.getFullPath(SpectrumPath)};
            Process p = Runtime.getRuntime().exec(msconvertcmd);
            LogManager.getRootLogger().info("MGF file coversion by msconvert.exe...." + SpectrumPath);
            LogManager.getRootLogger().debug("Command: " + Arrays.toString(msconvertcmd));
            PrintThread printThread = new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                LogManager.getRootLogger().info("msconvert : " + SpectrumPath + " failed");
                //PrintOutput(p);
                return;
            }
        } catch (IOException | InterruptedException ex) {
            java.util.logging.Logger.getLogger(MSConvert.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
}
