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

import java.io.BufferedReader;
import java.io.InputStreamReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class PrintThread extends Thread{

    Process process;
    
    public PrintThread(Process process){
        this.process=process;
    }
    
    @Override
    public void run() {
        try {
            final BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream()));
            String line = null;
            while ((line = reader.readLine()) != null) {                
                LogManager.getRootLogger().debug(line);
            }
            reader.close();
        } catch (final Exception e) {
            LogManager.getRootLogger().debug(e.getMessage());
        }        
    }
}
