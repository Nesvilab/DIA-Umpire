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
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.LCMSPeakStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.LogManager;

/**
 * For a given isolation window, extract grouped fragments for all peptide ions in the isolation window 
 * @author Chih-Chiang Tsou
 */
public class DIA_window_Quant implements  Runnable{

    public LCMSPeakDIAMS2 DIAWindow;
    HashMap<Integer, Integer> ScanClusterMap_Q1;
    HashMap<Integer, Integer> ScanClusterMap_Q2;
    HashMap<Integer, String> ScanClusterMap_Q3;
    String Q1Name;
    String Q2Name;    
    String Q3Name;
    LCMSPeakMS1 ms1lcms;
    LCMSID IDsummary;
    int NoThread=1;
    
    public DIA_window_Quant(String Q1Name,String Q2Name,String Q3Name,HashMap<Integer, Integer> ScanClusterMap_Q1, HashMap<Integer, Integer> ScanClusterMap_Q2, HashMap<Integer, String> ScanClusterMap_Q3, LCMSPeakMS1 ms1lcms,LCMSPeakDIAMS2 DIAWindow,LCMSID IDsummary, int NoThreads){
        this.ScanClusterMap_Q1=ScanClusterMap_Q1;
        this.ScanClusterMap_Q2=ScanClusterMap_Q2;
        this.ScanClusterMap_Q3=ScanClusterMap_Q3;
        this.Q1Name=Q1Name;
        this.Q2Name=Q2Name;
        this.Q3Name=Q3Name;
        this.ms1lcms=ms1lcms;
        this.DIAWindow=DIAWindow;
        this.IDsummary=IDsummary;
        this.NoThread=NoThreads;
    }
    @Override
    public void run() {
       
        if(!DIAWindow.ReadPeakCluster()){
            LogManager.getRootLogger().error("Reading Peak cluster result for " + DIAWindow.ScanCollectionName + " failed");
            return;
        }
        ExecutorService executorPool;
        executorPool = Executors.newFixedThreadPool(NoThread);
        
        //For each identified peptide ion, extract the precursor feature and grouped fragments from the isolation window
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            if (DIAWindow.DIA_MZ_Range.getX() <= pepIonID.GetPeakMz(2)&& DIAWindow.DIA_MZ_Range.getY() >= pepIonID.ObservedMz) {
                DIAMapClusterUnit mapunit = new DIAMapClusterUnit(pepIonID, Q1Name, Q2Name, Q3Name, ScanClusterMap_Q1, ScanClusterMap_Q2, ScanClusterMap_Q3, ms1lcms, DIAWindow);
                executorPool.execute(mapunit);
            }
        }
        executorPool.shutdown();

        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            LogManager.getRootLogger().info("interrupted..");
        }
        if (DIAWindow.datattype != SpectralDataType.DataType.pSMART) {
            if (!DIAWindow.ReadPrecursorFragmentClu2Cur()) {
                LogManager.getRootLogger().error("Reading precursor-fragment results for " + DIAWindow.ScanCollectionName + " failed");
                return;
            }

            for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
                for (PeakCluster cluster : pepIonID.MS1PeakClusters) {
                    if (DIAWindow.DIA_MZ_Range.getX() <= cluster.GetMaxMz() && DIAWindow.DIA_MZ_Range.getY() >= cluster.TargetMz()) {
                        DIAWindow.ExtractFragmentForPeakCluser(cluster);
                    }
                }
                for (PeakCluster ms2cluster : pepIonID.MS2UnfragPeakClusters) {
                    if (DIAWindow.DIA_MZ_Range.getX() <= ms2cluster.TargetMz() && DIAWindow.DIA_MZ_Range.getY() >= ms2cluster.TargetMz() && DIAWindow.PeakClusters.size()>=ms2cluster.Index) {
                        PeakCluster cluster = DIAWindow.PeakClusters.get(ms2cluster.Index - 1);
                        if (cluster.TargetMz() == ms2cluster.TargetMz() || cluster.Charge == ms2cluster.Charge) {
                            DIAWindow.ExtractFragmentForUnfragPeakCluser(cluster);
                        }
                    }
                }
            }
        }
        DIAWindow.ClearAllPeaks();
    }
    
}
