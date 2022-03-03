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
